"""Mapeamento INP <-> DOLFINx e montagem de cargas/reacoes no perfil FC.

So deve ser importado dentro da imagem com DOLFINx. Serial obrigatorio no FC:
mapas de nos, CLOAD e residual de RF nao foram validados com MPI/ghosts.
"""

from __future__ import annotations

import os
import time
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


def rss_mib() -> Tuple[float, float]:
    """(VmRSS, VmHWM) em MiB via /proc; (0, 0) se indisponivel."""
    rss_k = hwm_k = 0
    try:
        with open("/proc/self/status", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("VmRSS:"):
                    rss_k = int(line.split()[1])
                elif line.startswith("VmHWM:"):
                    hwm_k = int(line.split()[1])
    except OSError:
        return 0.0, 0.0
    return rss_k / 1024.0, hwm_k / 1024.0


def mem_profile_enabled() -> bool:
    return os.environ.get("FC_MEM_PROFILE", "").strip() not in ("", "0", "false", "False")


def log_fc_stage(name: str, t0: float, records: Optional[List[Dict[str, Any]]] = None) -> Dict[str, Any]:
    rss, hwm = rss_mib()
    rec = {
        "stage": name,
        "dt_s": round(time.perf_counter() - t0, 3),
        "rss_MiB": round(rss, 1),
        "hwm_MiB": round(hwm, 1),
    }
    if records is not None:
        records.append(rec)
    if mem_profile_enabled():
        print("FC_MEM {0}".format(rec), flush=True)
    return rec


_NEIGH27 = tuple((i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1))


class CoordLookup:
    """Indice O(N) de coordenadas com as mesmas tolerancias do mapa geometrico FC.

    Busca estrita em atol (1e-9) e, se vazio, folgada em loose (1e-6).
    Varios indices no raio sao ambiguos (coordenadas duplicadas).
    """

    def __init__(self, coords, atol: float = 1e-9, loose: float = 1e-6):
        self.coords = np.ascontiguousarray(coords, dtype=np.float64)
        if self.coords.ndim != 2 or self.coords.shape[1] < 3:
            raise ValueError("CoordLookup exige array (n, >=3), recebido {0}".format(self.coords.shape))
        self.coords = self.coords[:, :3]
        self.atol = float(atol)
        self.loose = float(loose)
        self._strict = self._build(self.atol)
        self._loose = self._build(self.loose)

    def _build(self, scale: float) -> Dict[Tuple[int, int, int], List[int]]:
        inv = 1.0 / scale
        keys = np.rint(self.coords * inv).astype(np.int64)
        table: Dict[Tuple[int, int, int], List[int]] = {}
        for i, key in enumerate(map(tuple, keys)):
            table.setdefault(key, []).append(i)
        return table

    def _hits(self, target: np.ndarray, radius: float, table) -> List[int]:
        key = tuple(int(v) for v in np.rint(target / radius))
        cand: List[int] = []
        for d in _NEIGH27:
            cand.extend(table.get((key[0] + d[0], key[1] + d[1], key[2] + d[2]), ()))
        if not cand:
            return []
        pts = self.coords[np.asarray(cand, dtype=np.int64)]
        dist = np.linalg.norm(pts - target, axis=1)
        return [cand[i] for i, d in enumerate(dist) if d <= radius]

    def find(self, xyz) -> Tuple[Optional[int], str, float]:
        """Retorna (indice, status, dist). status: ok | missing | ambiguous."""
        target = np.asarray(xyz, dtype=np.float64)[:3]
        hits = self._hits(target, self.atol, self._strict)
        radius = self.atol
        if not hits:
            hits = self._hits(target, self.loose, self._loose)
            radius = self.loose
        if not hits:
            # min dist so para a mensagem de erro (O(N), caminho raro)
            if self.coords.size:
                mind = float(np.linalg.norm(self.coords - target, axis=1).min())
            else:
                mind = float("nan")
            return None, "missing", mind
        if len(hits) > 1:
            pts = self.coords[np.asarray(hits, dtype=np.int64)]
            mind = float(np.linalg.norm(pts - target, axis=1).min())
            return None, "ambiguous", mind
        dist = float(np.linalg.norm(self.coords[hits[0]] - target))
        return int(hits[0]), "ok", dist


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
    lookup = CoordLookup(coords, atol=atol, loose=1e-6)
    mapping: Dict[int, int] = {}
    used: Dict[int, int] = {}
    for nid, node in model.nodes.items():
        target = (node.x, node.y, node.z)
        gi, status, mind = lookup.find(target)
        if status == "missing":
            raise FCProfileError(
                "No INP {0} {1} nao encontrado na geometria DOLFINx (min dist={2}).".format(
                    nid, [node.x, node.y, node.z], mind
                )
            )
        if status == "ambiguous":
            raise FCProfileError(
                "No INP {0} casa com varios vertices DOLFINx (coordenadas duplicadas); "
                "mapeamento ambiguo (min dist={1}).".format(nid, mind)
            )
        gi = int(gi)
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


class ComponentDofLocator:
    """Um collapse por componente e indice de coordenadas; lookup O(1) por no.

    Equivale a locate_dofs_geometrical por ponto, sem alocar FunctionSpace
    nem Function por no/BC/CLOAD.
    """

    def __init__(self, V, atol: float = 1e-9, loose: float = 1e-6):
        self.V = V
        self.gdim = int(V.mesh.geometry.dim)
        self.vcol = []
        self.parent_map = []
        self.lookup = []
        for comp in range(self.gdim):
            vcol, pmap = V.sub(comp).collapse()
            # dolfinx 0.11 devolve list[ndarray] com um bloco; precisa ser 1D
            # indexado pelo dof do subespaco (collapsed_dof -> parent_dof).
            pmap = np.asarray(pmap, dtype=np.int32).reshape(-1)
            coords = np.asarray(vcol.tabulate_dof_coordinates(), dtype=np.float64)
            if coords.ndim != 2:
                coords = np.reshape(coords, (-1, 3))
            if pmap.size != coords.shape[0]:
                raise RuntimeError(
                    "collapse map size={0} != n_dofs do subespaco {1}".format(
                        pmap.size, coords.shape[0]
                    )
                )
            self.vcol.append(vcol)
            self.parent_map.append(pmap)
            self.lookup.append(CoordLookup(coords, atol=atol, loose=loose))

    def locate(self, node, component: int):
        """DOFs do espaco pai V para um no INP e um componente 0..2.

        Retorno compativel com locate_component_parent_dofs: (parent, pair, Vcol).
        """
        if component < 0 or component >= self.gdim:
            return np.array([], dtype=np.int32), np.array([], dtype=np.int32), None
        idx, status, _mind = self.lookup[component].find((node.x, node.y, node.z))
        vcol = self.vcol[component]
        if status != "ok":
            return np.array([], dtype=np.int32), (np.array([], dtype=np.int32), np.array([], dtype=np.int32)), vcol
        sub = np.array([int(idx)], dtype=np.int32)
        parent = np.array([int(self.parent_map[component][int(idx)])], dtype=np.int32)
        return parent, (parent, sub), vcol


def component_dof_locator(V, atol: float = 1e-9, loose: float = 1e-6) -> ComponentDofLocator:
    loc = getattr(V, "_fc_comp_locator", None)
    if loc is None or getattr(loc, "V", None) is not V:
        loc = ComponentDofLocator(V, atol=atol, loose=loose)
        try:
            V._fc_comp_locator = loc
        except Exception:
            pass
    return loc


def locate_component_parent_dofs(V, node, component: int):
    """DOFs do espaco pai V para um no INP e um componente 0..2.

    Usa um ComponentDofLocator cacheado em V (um collapse por componente).
    """
    return component_dof_locator(V).locate(node, component)


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
    locator = component_dof_locator(V) if (step.loads or profile_fc) else None
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
        parent, _pair, _vcol = locator.locate(node, comp)
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
    """Lista (inp_node, component, parent_dof) apenas dos DOFs prescritos.

    Um Function + dirichletbc por (componente, valor), nunca por no. O valor
    prescrito que permanece e o ultimo *BOUNDARY daquele no/componente.
    """
    from dolfinx import default_scalar_type, fem

    rows = []
    seen_dof = set()
    seen_nc = set()
    last_value = {}
    last_pair = {}
    locator = component_dof_locator(V)
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
        parent, pair, Vcol = locator.locate(node, component)
        if parent.size == 0:
            if profile_fc:
                raise FCProfileError(
                    "BOUNDARY no INP {0} dof {1} nao localizou DOF no espaco DOLFINx.".format(
                        bc.node, bc.dof
                    )
                )
            continue
        key_nc = (int(bc.node), int(component))
        # ultimo valor prescrito vence no BC; o dof entra uma vez na tabela de reacoes
        last_value[key_nc] = float(bc.value)
        last_pair[key_nc] = (parent, pair, Vcol, component)
        if key_nc in seen_nc:
            continue
        seen_nc.add(key_nc)
        for dof in parent.tolist():
            if int(dof) in seen_dof and profile_fc:
                # dois nos distintos no mesmo dof nao deve ocorrer em serial
                pass
            seen_dof.add(int(dof))
            rows.append((int(bc.node), int(component), int(dof)))

    groups: Dict[Tuple[int, float], List[Tuple[np.ndarray, Any, Any]]] = {}
    for key_nc, value in last_value.items():
        parent, pair, Vcol, component = last_pair[key_nc]
        groups.setdefault((component, value), []).append((parent, pair, Vcol))

    bcs = []
    for (component, value), items in groups.items():
        parents = np.concatenate([it[0] for it in items]).astype(np.int32, copy=False)
        # pair = (parent_dofs, collapsed_dofs)
        sub_parts = []
        for _parent, pair, _vcol in items:
            if isinstance(pair, (tuple, list)) and len(pair) == 2:
                sub_parts.append(np.asarray(pair[1], dtype=np.int32))
            else:
                sub_parts.append(np.asarray(pair, dtype=np.int32))
        subs = np.concatenate(sub_parts).astype(np.int32, copy=False)
        Vcol = items[0][2]
        u_comp = fem.Function(Vcol)
        u_comp.x.array[:] = default_scalar_type(value)
        bcs.append(fem.dirichletbc(u_comp, (parents, subs), V.sub(component)))
    return bcs, rows


def residual_on_constrained(
    f_int, f_ext, constrained_rows, coords_by_nid: Dict[int, Tuple[float, float, float]]
) -> Dict[str, Any]:
    """R = f_int - f_ext somente nos DOFs prescritos.

    f_int deve vir da forma bilinear sem BC (assemble_vector de K u, esparso).
    Livres devem ter residual ~0. Nao monta K denso.
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
