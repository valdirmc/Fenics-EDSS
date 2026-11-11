#!/usr/bin/env python3
"""
Executa análise a partir de um arquivo .inp.
Detecta o tipo (estática ou térmica) e chama a função apropriada.
"""

import sys
import os
import numpy as np
from mpi4py import MPI
from dolfinx import mesh, fem, default_scalar_type
from dolfinx.fem.petsc import assemble_matrix, assemble_vector, apply_lifting, set_bc
from petsc4py import PETSc
import ufl
import basix.ufl

sys.path.insert(0, "/app/src")
from inp_parser import parse_inp, summary, Step, select_analysis, get_material_for_element
from job_io import (
    USAGE,
    parse_run_args,
    prepare_output_dir,
    resolve_output_dir,
    step_output_base,
    write_run_manifest,
    write_vtk_function,
    verify_run_artifacts,
)
from fc_profile import (
    FCProfileError,
    check_honored_outputs,
    normalize_units,
    required_output_fields,
)


def build_mesh_from_inp(model, debug=0):
    """Constrói uma malha DOLFINx a partir do modelo INP."""
    if not model.elements:
        raise ValueError("Nenhum elemento encontrado no .inp")

    first_el = next(iter(model.elements.values()))
    etype = first_el.type.upper()

    node_ids = sorted(model.nodes.keys())
    node_map = {nid: i for i, nid in enumerate(node_ids)}

    coords = np.array([[model.nodes[nid].x, model.nodes[nid].y, model.nodes[nid].z]
                       for nid in node_ids], dtype=np.float64)

    cells = []
    for el in model.elements.values():
        if el.type.upper() != etype:
            raise ValueError(f"Elementos mistos não suportados ({el.type} vs {etype})")
        nodes = [node_map[n] for n in el.nodes]
        if etype == "C3D8" and len(nodes) == 8:
            nodes = [nodes[0], nodes[1], nodes[3], nodes[2],
                     nodes[4], nodes[5], nodes[7], nodes[6]]
        cells.append(nodes)
    cells = np.array(cells, dtype=np.int32)

    if etype == "C3D8":
        cell_name = "hexahedron"
    elif etype == "C3D4":
        cell_name = "tetrahedron"
    else:
        raise ValueError(f"Tipo de elemento não suportado: {etype}")

    if debug:
        print(f"Malha: {len(node_ids)} nós, {len(cells)} elementos ({etype})")

    coord_element = basix.ufl.element("Lagrange", cell_name, 1, shape=(coords.shape[1],))
    ufl_mesh = ufl.Mesh(coord_element)
    domain = mesh.create_mesh(MPI.COMM_WORLD, cells, ufl_mesh, coords)
    return domain, node_map


def run_linear_static_step(model, step, output_base: str, debug: int = 0, profile: str = None):
    """Executa um step de estatica linear com multiplos materiais."""
    profile_fc = (profile or "").strip().lower() == "fc"
    domain, node_map = build_mesh_from_inp(model, debug=debug)
    
    domain.topology.create_connectivity(0, domain.topology.dim)
    domain.topology.create_connectivity(domain.topology.dim, 0)
    
    V = fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim,)))

    # Materiais por célula (DG-0)
    Q = fem.functionspace(domain, ("DG", 0))
    E_func = fem.Function(Q)
    nu_func = fem.Function(Q)
    
    # Inicializa com o primeiro material
    if not model.materials:
        raise ValueError("Nenhum material *ELASTIC encontrado no .inp")
    default_mat = next(iter(model.materials.values()))
    if profile_fc and (default_mat.E <= 0):
        raise FCProfileError("Perfil FC exige *ELASTIC com E e nu informados (sem default silencioso).")
    E_func.x.array[:] = default_mat.E if default_mat.E > 0 else 210000
    nu_func.x.array[:] = default_mat.nu if default_mat.nu > 0 else 0.3
    
    # Sobrescreve com materiais específicos por elset
    elem_to_cell = {el.id: i for i, el in enumerate(model.elements.values())}
    
    for elset_name, elem_ids in model.elsets.items():
        mat = get_material_for_element(model, elem_ids[0] if elem_ids else 0)
        if mat is None:
            continue
        for eid in elem_ids:
            if eid in elem_to_cell:
                cell_idx = elem_to_cell[eid]
                if cell_idx < len(E_func.x.array):
                    if profile_fc and mat.E <= 0:
                        raise FCProfileError("Material sem E no perfil FC.")
                    E_func.x.array[cell_idx] = mat.E if mat.E > 0 else 210000
                    nu_func.x.array[cell_idx] = mat.nu if mat.nu > 0 else 0.3
    
    E_func.x.scatter_forward()
    nu_func.x.scatter_forward()
    
    if debug:
        print(f"Materiais aplicados: E min={E_func.x.array.min():.0f}, max={E_func.x.array.max():.0f}")

    # Propriedades como funções
    E = E_func
    nu = nu_func
    mu = E / (2.0 * (1.0 + nu))
    lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))

    # BCs por componente, honrando o valor prescrito (nao zerar tudo).
    bcs = []
    for bc in model.boundary:
        if bc.node not in model.nodes:
            continue
        n = model.nodes[bc.node]
        component = bc.dof - 1
        if component < 0 or component >= domain.geometry.dim:
            continue

        def marker(x, nx=n.x, ny=n.y, nz=n.z):
            return (np.isclose(x[0], nx, atol=1e-6) &
                    np.isclose(x[1], ny, atol=1e-6) &
                    np.isclose(x[2], nz, atol=1e-6))

        Vcol, _map = V.sub(component).collapse()
        dofs = fem.locate_dofs_geometrical((V.sub(component), Vcol), marker)
        if len(dofs) > 0:
            u_comp = fem.Function(Vcol)
            u_comp.x.array[:] = bc.value
            bcs.append(fem.dirichletbc(u_comp, dofs, V.sub(component)))

    # Formulação
    def epsilon(u):
        return ufl.sym(ufl.grad(u))

    def sigma(u):
        return lmbda * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2.0 * mu * epsilon(u)

    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx
    L_form = ufl.dot(fem.Constant(domain, default_scalar_type((0.0, 0.0, 0.0))), v) * ufl.dx

    a_form = fem.form(a)
    L_compiled = fem.form(L_form)
    A = assemble_matrix(a_form, bcs=bcs)
    A.assemble()
    b = assemble_vector(L_compiled)

    # Cargas do step
    for cl in step.loads:
        if cl.node not in model.nodes:
            continue
        n = model.nodes[cl.node]
        dist = np.linalg.norm(coords - np.array([n.x, n.y, n.z]), axis=1)
        node_idx = np.argmin(dist)
        if dist[node_idx] > 1e-6:
            continue
        comp = cl.dof - 1
        dof_comp = node_idx * 3 + comp
        with b.localForm() as local:
            if dof_comp < local.size:
                local[dof_comp] += cl.value

    # DLOAD do step
    for dl in step.dloads:
        for elset_name, elem_ids in model.elsets.items():
            if elset_name.upper() == dl.elset.upper():
                for eid in elem_ids:
                    if eid in model.elements:
                        for nid in model.elements[eid].nodes:
                            if nid in model.nodes:
                                n = model.nodes[nid]
                                dist = np.linalg.norm(coords - np.array([n.x, n.y, n.z]), axis=1)
                                node_idx = np.argmin(dist)
                                if dist[node_idx] < 1e-6:
                                    dof_z = node_idx * 3 + 2
                                    with b.localForm() as local:
                                        if dof_z < local.size:
                                            local[dof_z] += dl.value * 0.05
                break

    f_ext = b.copy()

    apply_lifting(b, [a_form], [bcs])
    b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
    set_bc(b, bcs)

    solver = PETSc.KSP().create(domain.comm)
    solver.setOperators(A)
    solver.setType("preonly")
    solver.getPC().setType("lu")
    uh = fem.Function(V)
    solver.solve(b, uh.x.petsc_vec)
    uh.x.scatter_forward()
    uh.name = "displacement"

    write_vtk_function(domain.comm, uh, output_base, debug=debug)
    honored = ["U"]
    artifacts = [{"path": output_base + ".pvd", "field": "U", "location": "nodal-CG1"}]
    extra = {"honored": honored, "artifacts": artifacts, "reactions": None}

    if profile_fc:
        extra = _write_fc_linear_fields(
            domain, uh, v, sigma, epsilon, f_ext, a_form, model, output_base, debug
        )

    norma = np.linalg.norm(uh.x.array)
    if debug:
        print(f"Norma: {norma:.6e}")
    return uh, extra

def run_from_inp(inp_path: str, output_dir: str = "/app/output", debug: int = 0,
                 profile: str = None, units: str = None):
    """Executa cada step e grava a saida somente em output_dir."""
    if not os.path.exists(inp_path):
        print(f"ERRO: arquivo nao encontrado: {inp_path}")
        sys.exit(1)

    profile_fc = (profile or "").strip().lower() == "fc"
    units_norm = None
    if profile_fc:
        if not units:
            raise FCProfileError("Perfil FC exige --units (use mm-N-MPa-t-s).")
        units_norm = normalize_units(units)

    output_dir = prepare_output_dir(resolve_output_dir(output_dir))
    model = parse_inp(inp_path, debug=debug, profile=profile)
    if debug:
        print(summary(model))
        print("---")

    if not model.steps:
        model.steps.append(Step())

    results = []
    for i, step in enumerate(model.steps):
        step_num = i + 1
        output_base = step_output_base(output_dir, inp_path, step_num)

        if debug:
            print(f"\n=== Step {step_num} ===")

        kind = select_analysis(model, step)
        if profile_fc and kind != "linear":
            raise FCProfileError(
                "Perfil FC aceita somente estatica linear (recebido: {0}).".format(kind)
            )
        if kind == "heat":
            from run_thermal import run_thermal
            result = run_thermal(inp_path, output_base, debug)
            extra = {"honored": [], "artifacts": [], "reactions": None}
        elif kind == "hyper":
            from run_hyper import run_hyperelastic
            result = run_hyperelastic(model, step, output_base, debug)
            extra = {"honored": [], "artifacts": [], "reactions": None}
        else:
            result, extra = run_linear_static_step(
                model, step, output_base, debug, profile=profile
            )
        results.append((result, extra))

    if profile_fc:
        _finalize_fc_run(model, output_dir, inp_path, units_norm, results, debug)

    return results




def _write_fc_linear_fields(domain, uh, v, sigma, epsilon, f_ext, a_form, model, output_base, debug):
    """Grava U, S, principais, von Mises, E e RF do perfil FC. Falha se VTK falhar."""
    import json
    import os
    from dolfinx import fem
    from fc_fields import (
        interpolate_cell_scalar,
        interpolate_cell_tensor,
        interpolate_cell_vector,
        principals_algebraic,
        tensor_array_from_function,
        write_functions_vtk,
    )

    uh.name = "displacement"
    strain_ufl = epsilon(uh)
    cauchy_ufl = sigma(uh)
    cauchy_f = interpolate_cell_tensor(domain, cauchy_ufl, "cauchy_stress")
    strain_f = interpolate_cell_tensor(domain, strain_ufl, "strain")
    sigma_cells = tensor_array_from_function(cauchy_f)
    principals, vm = principals_algebraic(sigma_cells)
    prin_f = interpolate_cell_vector(domain, principals, "principal_stress")
    vm_f = interpolate_cell_scalar(domain, vm, "von_mises")

    cell_pvd = write_functions_vtk(
        domain.comm,
        [cauchy_f, strain_f, prin_f, vm_f],
        output_base + "_cell",
    )

    # Reacoes: R = K u - f_ext nos dofs (nos presos devem equilibrar as cargas).
    f_int = assemble_vector(fem.form(ufl.inner(sigma(uh), epsilon(v)) * ufl.dx))
    try:
        f_int.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
    except Exception:
        pass
    rvec = f_int.array.copy()
    try:
        ext = f_ext.array
        n = min(rvec.size, ext.size)
        rvec[:n] = rvec[:n] - ext[:n]
    except Exception:
        # petsc4py Vec
        ext = f_ext.getArray(readonly=True)
        n = min(rvec.size, ext.size)
        rvec[:n] = rvec[:n] - ext[:n]

    V = uh.function_space
    rf = fem.Function(V)
    rf.name = "reaction_force"
    rf.x.array[:] = 0.0
    mlen = min(rf.x.array.size, rvec.size)
    rf.x.array[:mlen] = rvec[:mlen]
    rf.x.scatter_forward()
    rf_pvd = write_functions_vtk(domain.comm, [rf], output_base + "_RF")

    totals = {"RFx": 0.0, "RFy": 0.0, "RFz": 0.0}
    # Soma nos nos de model.boundary (presos)
    coords = domain.geometry.x
    dim = domain.geometry.dim
    seen = set()
    for bc in model.boundary:
        if bc.node in seen or bc.node not in model.nodes:
            continue
        seen.add(bc.node)
        n = model.nodes[bc.node]
        dist = np.linalg.norm(coords - np.array([n.x, n.y, n.z]), axis=1)
        node_idx = int(np.argmin(dist))
        if dist[node_idx] > 1e-6:
            continue
        totals["RFx"] += float(rvec[node_idx * dim + 0]) if node_idx * dim + 0 < rvec.size else 0.0
        totals["RFy"] += float(rvec[node_idx * dim + 1]) if node_idx * dim + 1 < rvec.size else 0.0
        totals["RFz"] += float(rvec[node_idx * dim + 2]) if node_idx * dim + 2 < rvec.size else 0.0

    react_path = output_base + "_reactions.json"
    with open(react_path, "w", encoding="utf-8") as handle:
        json.dump({"totals": totals, "n_constrained_nodes": len(seen)}, handle, indent=2)
        handle.write("\n")

    if debug:
        print("FC campos celula:", cell_pvd)
        print("FC RF:", totals)

    honored = ["U", "S", "E", "PRINCIPAL", "VON_MISES", "RF"]
    artifacts = [
        {"path": output_base + ".pvd", "field": "U", "location": "nodal-CG1", "space": "Lagrange-P1"},
        {"path": cell_pvd, "field": "S,E,PRINCIPAL,VON_MISES", "location": "cell-DG0-constant-C3D4-P1",
         "tensor_order": "row-major XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ",
         "principals": "algebraic sigma1>=sigma2>=sigma3 with sign",
         "note": "nodal smoothed does not replace cell field"},
        {"path": rf_pvd, "field": "RF", "location": "nodal"},
        {"path": react_path, "field": "RF_totals", "location": "json"},
    ]
    return {"honored": honored, "artifacts": artifacts, "reactions": totals}


def _stack_versions():
    info = {"python": sys.version.split()[0]}
    for name in ("dolfinx", "ufl", "basix", "ffcx", "petsc4py", "mpi4py", "numpy"):
        try:
            mod = __import__(name)
            info[name] = getattr(mod, "__version__", "?")
        except Exception:
            info[name] = "ausente"
    return info


def _finalize_fc_run(model, output_dir, inp_path, units_norm, results, debug):
    rec = model.fc_record
    honored = []
    artifacts = []
    reactions = None
    for _res, extra in results:
        honored = list(extra.get("honored") or [])
        artifacts.extend(extra.get("artifacts") or [])
        reactions = extra.get("reactions") or reactions

    requested = required_output_fields(rec) if rec else ["U", "S", "PRINCIPAL", "VON_MISES"]
    check_honored_outputs(requested, honored)

    payload = {
        "profile": "fc",
        "units": units_norm,
        "inp": os.path.abspath(inp_path),
        "start_time": 0.0,
        "n_steps": len(model.steps),
        "element_types": sorted({el.type.upper() for el in model.elements.values()}),
        "n_nodes": len(model.nodes),
        "n_elements": len(model.elements),
        "materials": {
            name: {"E": mat.E, "nu": mat.nu} for name, mat in model.materials.items()
        },
        "ignored_cards": [
            {"keyword": it.keyword, "reason": it.reason, "line": it.line, "detail": it.detail}
            for it in (rec.ignored if rec else [])
        ],
        "output_requests": [
            {"kind": r.kind, "fields": r.fields, "nset": r.nset, "totals_only": r.totals_only, "line": r.line}
            for r in (rec.requests if rec else [])
        ],
        "honored_fields": honored,
        "field_contract": {
            "U": "nodal Lagrange P1",
            "S": "Cauchy cell DG0, constant on C3D4 P1; not nodal-smoothed",
            "E": "infinitesimal strain cell DG0",
            "PRINCIPAL": "algebraic sigma1>=sigma2>=sigma3 with sign, cell DG0",
            "VON_MISES": "cell DG0",
            "RF": "reactions at constrained nodes and totals",
        },
        "static_increments": rec.static_increments if rec else None,
        "reactions": reactions,
        "solver_versions": _stack_versions(),
        "image": {
            "name_requested": os.environ.get("FENICS_IMAGE", "edss/fenicsx:dev"),
            "note": (
                "Versoes acima sao as importadas nesta execucao. "
                "Image ID local (docker image inspect .Id) identifica o artefato nesta maquina "
                "e NAO e um RepoDigest verificavel de registry, nem referencia de pull. "
                "Nao usar <repo>:<tag>@sha256:<imageId> como pin de registry."
            ),
        },
        "artifacts": artifacts,
        "limits": {
            "wheel_validation": False,
            "independent_benchmark": False,
            "material_plan_70000": "nao corrigido automaticamente se o .inp trouxer outro E",
        },
    }
    manifest_path = write_run_manifest(output_dir, payload)
    required = [os.path.basename(a["path"]) for a in artifacts if a.get("path")]
    required.append("run_manifest.json")
    # include manifest itself in verification after writing hash
    # re-write payload with manifest listed then verify others
    verify_run_artifacts(manifest_path, required[:-1])
    if debug:
        print("Manifesto:", manifest_path)


if __name__ == "__main__":
    try:
        args = parse_run_args(sys.argv[1:])
    except ValueError as exc:
        print(f"ERRO: {exc}")
        print(USAGE)
        sys.exit(2)
    if not args.inp_file:
        print(USAGE)
        sys.exit(1)
    try:
        run_from_inp(
            args.inp_file,
            output_dir=resolve_output_dir(args.output_dir),
            debug=args.debug,
            profile=args.profile,
            units=args.units,
        )
    except FCProfileError as exc:
        print("ERRO FC: {0}".format(exc))
        sys.exit(2)
    except RuntimeError as exc:
        print("ERRO: {0}".format(exc))
        sys.exit(3)
    except Exception as exc:
        print("ERRO: {0}: {1}".format(type(exc).__name__, exc))
        sys.exit(3)
