"""Mapeamento INP <-> DOLFINx e montagem de cargas/reacoes no perfil FC.

So deve ser importado dentro da imagem com DOLFINx. Serial obrigatorio no FC:
mapas de nos, CLOAD e residual de RF nao foram validados com MPI/ghosts.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from fc_profile import FCProfileError


def require_serial_fc(comm, profile_fc: bool) -> None:
    if not profile_fc:
        return
    size = int(comm.size)
    if size > 1:
        raise FCProfileError(
            "Perfil FC neste lote so e suportado em serial (MPI size=1). "
            "size={0}: indexacao de nos/CLOAD/RF com ghosts nao foi validada.".format(size)
        )


def _geom_dofmap_links(gdmap, cell: int) -> List[int]:
    try:
        return [int(i) for i in gdmap[cell]]
    except Exception:
        return [int(i) for i in gdmap.links(cell)]


def inp_node_to_geom(domain, model, atol: float = 1e-9) -> Dict[int, int]:
    """INP node id -> indice de vertice em domain.geometry.x (apos create_mesh).

    create_mesh pode reordenar vertices; o node_map do array de entrada
    nao e o indice de DOF. Coordenadas duplicadas ou no ausente sao erro.
    """
    coords = np.asarray(domain.geometry.x)
    mapping: Dict[int, int] = {}
    used: Dict[int, int] = {}
    for nid, node in model.nodes.items():
        target = np.array([node.x, node.y, node.z], dtype=float)
        dist = np.linalg.norm(coords - target, axis=1)
        hits = np.flatnonzero(dist <= atol)
        if hits.size == 0:
            hits = np.flatnonzero(dist <= 1e-6)
        if hits.size == 0:
            raise FCProfileError(
                "No INP {0} {1} nao encontrado na geometria DOLFINx (min dist={2}).".format(
                    nid, target.tolist(), float(dist.min()) if dist.size else "n/a"
                )
            )
        if hits.size > 1:
            raise FCProfileError(
                "No INP {0} casa com {1} vertices DOLFINx (coordenadas duplicadas); "
                "mapeamento ambiguo.".format(nid, int(hits.size))
            )
        gi = int(hits[0])
        if gi in used and used[gi] != nid:
            raise FCProfileError(
                "Vertice DOLFINx {0} mapeado a dois nos INP ({1} e {2}).".format(
                    gi, used[gi], nid
                )
            )
        used[gi] = nid
        mapping[int(nid)] = gi
    return mapping


def inp_element_to_cell(domain, model, nid_to_geom: Dict[int, int]) -> Dict[int, int]:
    """INP element id -> indice de celula local via conjunto de vertices."""
    geom = domain.geometry
    if hasattr(geom, "dofmaps"):
        gdmap = geom.dofmaps[0]
    else:
        gdmap = geom.dofmap
    dim = domain.topology.dim
    ncells = int(domain.topology.index_map(dim).size_local)
    key_to_cell: Dict[Tuple[int, ...], int] = {}
    for c in range(ncells):
        key = tuple(sorted(_geom_dofmap_links(gdmap, c)))
        if key in key_to_cell:
            raise FCProfileError(
                "Duas celulas DOLFINx com os mesmos vertices: {0}".format(key)
            )
        key_to_cell[key] = c
    mapping: Dict[int, int] = {}
    for eid, el in model.elements.items():
        try:
            key = tuple(sorted(nid_to_geom[n] for n in el.nodes))
        except KeyError as exc:
            raise FCProfileError(
                "Elemento INP {0} referencia no sem mapa geometrico: {1}".format(eid, exc)
            ) from exc
        if key not in key_to_cell:
            raise FCProfileError(
                "Elemento INP {0} (vertices geom {1}) nao encontrado nas celulas DOLFINx.".format(
                    eid, key
                )
            )
        mapping[int(eid)] = key_to_cell[key]
    if len(set(mapping.values())) != len(mapping):
        raise FCProfileError("Dois elementos INP cairam na mesma celula DOLFINx.")
    if len(mapping) != ncells:
        raise FCProfileError(
            "Mapa INP->celula incompleto ({0} elementos, {1} celulas).".format(
                len(mapping), ncells
            )
        )
    return mapping


def locate_component_parent_dofs(V, node, component: int):
    """DOFs do espaco pai V para um no INP e um componente 0..2."""
    from dolfinx import fem

    def marker(x, nx=node.x, ny=node.y, nz=node.z):
        return (
            np.isclose(x[0], nx, atol=1e-9)
            & np.isclose(x[1], ny, atol=1e-9)
            & np.isclose(x[2], nz, atol=1e-9)
        )

    Vcol, _map = V.sub(component).collapse()
    dofs = fem.locate_dofs_geometrical((V.sub(component), Vcol), marker)
    if isinstance(dofs, (tuple, list)) and len(dofs) == 2:
        parent = np.asarray(dofs[0], dtype=np.int32)
        pair = dofs
    else:
        parent = np.asarray(dofs, dtype=np.int32)
        pair = dofs
    return parent, pair, Vcol


def add_to_vector(vec, dof: int, value: float) -> None:
    dof = int(dof)
    try:
        with vec.localForm() as local:
            if dof < 0 or dof >= local.size:
                raise IndexError("dof {0} fora de localForm size={1}".format(dof, local.size))
            local[dof] += value
            return
    except Exception as first:
        if hasattr(vec, "array") and 0 <= dof < vec.array.size:
            vec.array[dof] += value
            return
        raise RuntimeError(
            "nao foi possivel acumular carga no dof {0}: {1}".format(dof, first)
        ) from first


def vec_as_array(vec) -> np.ndarray:
    if hasattr(vec, "array"):
        return np.asarray(vec.array).copy()
    if hasattr(vec, "getArray"):
        return np.asarray(vec.getArray()).copy()
    raise RuntimeError("vetor sem interface .array/getArray")


def apply_cloads(V, b, model, step, profile_fc: bool) -> Dict[str, float]:
    """Acumula CLOAD nos DOFs localizados por coordenada do no INP (nao node*3)."""
    applied = {"Fx": 0.0, "Fy": 0.0, "Fz": 0.0}
    keys = ("Fx", "Fy", "Fz")
    for cl in step.loads:
        if cl.node not in model.nodes:
            if profile_fc:
                raise FCProfileError("CLOAD referencia no inexistente: {0}".format(cl.node))
            continue
        node = model.nodes[cl.node]
        comp = cl.dof - 1
        if comp < 0 or comp >= V.mesh.geometry.dim:
            if profile_fc:
                raise FCProfileError("CLOAD com DOF invalido no no {0}: {1}".format(cl.node, cl.dof))
            continue
        parent, _pair, _vcol = locate_component_parent_dofs(V, node, comp)
        if parent.size == 0:
            if profile_fc:
                raise FCProfileError(
                    "CLOAD no INP {0} dof {1} nao localizou DOF no espaco DOLFINx.".format(
                        cl.node, cl.dof
                    )
                )
            continue
        if profile_fc and parent.size != 1:
            raise FCProfileError(
                "CLOAD no INP {0} dof {1} localizou {2} DOFs (esperado 1 em serial).".format(
                    cl.node, cl.dof, int(parent.size)
                )
            )
        for dof in parent.tolist():
            add_to_vector(b, int(dof), float(cl.value))
        applied[keys[comp]] += float(cl.value)
    return applied


def collect_constrained_dofs(V, model, profile_fc: bool):
    """Lista (inp_node, component, parent_dof) apenas dos DOFs prescritos."""
    from dolfinx import fem

    rows = []
    bcs = []
    seen_dof = set()
    for bc in model.boundary:
        if bc.node not in model.nodes:
            if profile_fc:
                raise FCProfileError("BOUNDARY referencia no inexistente: {0}".format(bc.node))
            continue
        node = model.nodes[bc.node]
        component = bc.dof - 1
        if component < 0 or component >= V.mesh.geometry.dim:
            if profile_fc:
                raise FCProfileError("BOUNDARY com DOF invalido no no {0}: {1}".format(bc.node, bc.dof))
            continue
        parent, pair, Vcol = locate_component_parent_dofs(V, node, component)
        if parent.size == 0:
            if profile_fc:
                raise FCProfileError(
                    "BOUNDARY no INP {0} dof {1} nao localizou DOF no espaco DOLFINx.".format(
                        bc.node, bc.dof
                    )
                )
            continue
        u_comp = fem.Function(Vcol)
        u_comp.x.array[:] = bc.value
        bcs.append(fem.dirichletbc(u_comp, pair, V.sub(component)))
        for dof in parent.tolist():
            key = (int(bc.node), int(component), int(dof))
            if (int(bc.node), int(component)) in {(r[0], r[1]) for r in rows}:
                # mesmo no/componente repetido: permanece o ultimo valor prescrito no BC,
                # mas o dof entra uma vez so na tabela de reacoes.
                continue
            if int(dof) in seen_dof and profile_fc:
                # dois nos distintos no mesmo dof nao deve ocorrer em serial
                pass
            seen_dof.add(int(dof))
            rows.append(key)
    return bcs, rows


def residual_on_constrained(
    f_int, f_ext, constrained_rows, coords_by_nid: Dict[int, Tuple[float, float, float]]
) -> Dict[str, Any]:
    """R = f_int - f_ext somente nos DOFs prescritos.

    f_int deve vir da forma bilinear sem BC (K_full u). Livres devem ter residual ~0.
    """
    r_all = vec_as_array(f_int)
    ext = vec_as_array(f_ext)
    n = min(r_all.size, ext.size)
    residual = np.zeros_like(r_all)
    residual[:n] = r_all[:n] - ext[:n]

    totals = {"RFx": 0.0, "RFy": 0.0, "RFz": 0.0}
    keys = ("RFx", "RFy", "RFz")
    moments = {"Mx": 0.0, "My": 0.0, "Mz": 0.0}
    rf_nonzero = np.zeros_like(residual)
    listed = []
    used = set()
    for nid, comp, dof in constrained_rows:
        if dof < 0 or dof >= residual.size:
            raise RuntimeError("DOF prescrito {0} fora do residual (size={1}).".format(dof, residual.size))
        if (nid, comp) in used:
            continue
        used.add((nid, comp))
        val = float(residual[dof])
        rf_nonzero[dof] = val
        totals[keys[comp]] += val
        x, y, z = coords_by_nid[nid]
        fx = val if comp == 0 else 0.0
        fy = val if comp == 1 else 0.0
        fz = val if comp == 2 else 0.0
        moments["Mx"] += y * fz - z * fy
        moments["My"] += z * fx - x * fz
        moments["Mz"] += x * fy - y * fx
        listed.append({"inp_node": int(nid), "component": int(comp), "dof": int(dof), "R": val})

    constrained_set = {dof for _n, _c, dof in constrained_rows}
    free_mask = np.ones(residual.size, dtype=bool)
    for dof in constrained_set:
        if 0 <= dof < residual.size:
            free_mask[dof] = False
    free_abs = np.abs(residual[free_mask]) if residual.size else np.array([0.0])
    return {
        "residual_all": residual,
        "rf_constrained": rf_nonzero,
        "totals": totals,
        "moment_R_origin": moments,
        "n_constrained_dofs": len(used),
        "residual_free_dofs_max_abs": float(free_abs.max()) if free_abs.size else 0.0,
        "entries": listed,
    }


def check_solution_finite(uh, solver=None) -> None:
    arr = np.asarray(uh.x.array)
    if arr.size == 0 or not np.isfinite(arr).all():
        raise RuntimeError("Solucao contem NaN/Inf ou esta vazia; recusado como sucesso.")
    if solver is not None:
        try:
            reason = int(solver.getConvergedReason())
        except Exception:
            return
        if reason <= 0:
            raise RuntimeError("PETSc KSP nao convergiu (reason={0}).".format(reason))
