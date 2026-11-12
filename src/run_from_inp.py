#!/usr/bin/env python3
"""
Executa análise a partir de um arquivo .inp.
Detecta o tipo (estática ou térmica) e chama a função apropriada.
"""

import json
import time
import sys
import os
import uuid
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
    write_run_status,
    write_vtk_function,
    verify_run_artifacts,
)
from fc_profile import (
    FCProfileError,
    check_honored_outputs,
    normalize_units,
    parse_report_payload,
    required_output_fields,
    PARSE_NOT_SOLVED_BANNER,
)
from fc_assembly import (
    apply_cloads,
    check_solution_finite,
    collect_constrained_dofs,
    inp_element_to_cell,
    inp_node_to_geom,
    log_fc_stage,
    mem_profile_enabled,
    require_serial_fc,
    residual_on_constrained,
)
from fc_solver import SOLVER_LU, solve_linear
from job_io import normalize_rtol, normalize_solver


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


def run_linear_static_step(model, step, output_base: str, debug: int = 0, profile: str = None,
                           solver: str = "lu", rtol: float = 1e-10):
    """Executa um step de estatica linear com multiplos materiais."""
    profile_fc = (profile or "").strip().lower() == "fc"
    solver = normalize_solver(solver)
    rtol = normalize_rtol(rtol, solver)
    require_serial_fc(MPI.COMM_WORLD, profile_fc)
    t0 = time.perf_counter()
    stages = [] if mem_profile_enabled() else None
    domain, node_map = build_mesh_from_inp(model, debug=debug)
    log_fc_stage("mesh", t0, stages)

    domain.topology.create_connectivity(0, domain.topology.dim)
    domain.topology.create_connectivity(domain.topology.dim, 0)

    nid_to_geom = inp_node_to_geom(domain, model)
    log_fc_stage("inp_node_to_geom", t0, stages)
    eid_to_cell = inp_element_to_cell(domain, model, nid_to_geom)
    log_fc_stage("inp_element_to_cell", t0, stages)
    if debug:
        print(f"Mapa INP->geom: {len(nid_to_geom)} nos; INP->celula: {len(eid_to_cell)} elems")

    V = fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim,)))
    log_fc_stage("functionspace", t0, stages)

    # Materiais por célula (DG-0). Indices de celula via mapa INP, nao enumeracao crua.
    Q = fem.functionspace(domain, ("DG", 0))
    E_func = fem.Function(Q)
    nu_func = fem.Function(Q)

    if not model.materials:
        raise ValueError("Nenhum material *ELASTIC encontrado no .inp")
    default_mat = next(iter(model.materials.values()))
    if profile_fc and (default_mat.E <= 0):
        raise FCProfileError("Perfil FC exige *ELASTIC com E e nu informados (sem default silencioso).")
    E_func.x.array[:] = default_mat.E if default_mat.E > 0 else 210000
    nu_func.x.array[:] = default_mat.nu if default_mat.nu > 0 else 0.3

    for eid in model.elements:
        mat = get_material_for_element(model, eid)
        if mat is None:
            continue
        cell_idx = eid_to_cell[eid]
        if cell_idx < len(E_func.x.array):
            if profile_fc and mat.E <= 0:
                raise FCProfileError("Material sem E no perfil FC.")
            E_func.x.array[cell_idx] = mat.E if mat.E > 0 else 210000
            nu_func.x.array[cell_idx] = mat.nu if mat.nu > 0 else 0.3

    E_func.x.scatter_forward()
    nu_func.x.scatter_forward()
    log_fc_stage("materials", t0, stages)

    if debug:
        print(f"Materiais aplicados: E min={E_func.x.array.min():.0f}, max={E_func.x.array.max():.0f}")

    E = E_func
    nu = nu_func
    mu = E / (2.0 * (1.0 + nu))
    lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))

    bcs, constrained_rows = collect_constrained_dofs(V, model, profile_fc)
    log_fc_stage("collect_constrained_dofs", t0, stages)
    if profile_fc and not constrained_rows:
        raise RuntimeError(
            "Perfil FC recusa sistema sem DOFs prescritos (*BOUNDARY). "
            "Ha modos de corpo rigido; uma fatoracao LU que nao acusa "
            "singularidade nao conta como sucesso."
        )

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

    # RF usa assemble_vector da forma (K u esparso), nao uma copia densa nem A_full.
    # Matriz com BC so no solve.
    A = assemble_matrix(a_form, bcs=bcs)
    A.assemble()
    log_fc_stage("assemble_A", t0, stages)
    b = assemble_vector(L_compiled)
    log_fc_stage("assemble_b", t0, stages)

    applied = apply_cloads(V, b, model, step, profile_fc)
    log_fc_stage("apply_cloads", t0, stages)

    if step.dloads:
        if profile_fc:
            raise FCProfileError("Perfil FC rejeita *DLOAD (sem heuristicas de pressao).")
        # Recorte antigo (nao FC): contribuicao nodal aproximada; nao usar como pressao real.
        coords = domain.geometry.x
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

    uh = fem.Function(V)
    solver_info = solve_linear(
        A, b, uh, V, solver=solver, rtol=rtol, debug=debug,
    )
    log_fc_stage("solve", t0, stages)
    uh.name = "displacement"

    write_vtk_function(domain.comm, uh, output_base, debug=debug)
    honored = ["U"]
    artifacts = [{"path": output_base + ".pvd", "field": "U", "location": "nodal-CG1"}]
    extra = {
        "honored": honored, "artifacts": artifacts, "reactions": None,
        "applied": applied, "solver_info": solver_info,
    }

    if profile_fc:
        fields = _write_fc_linear_fields(
            domain, uh, v, sigma, epsilon, f_ext, a_form, model, output_base, debug,
            constrained_rows=constrained_rows, applied=applied,
        )
        fields["solver_info"] = solver_info
        extra = fields
        log_fc_stage("fc_fields_rf", t0, stages)
        _check_fc_force_balance(extra)

    norma = np.linalg.norm(uh.x.array)
    if debug:
        print(f"Norma: {norma:.6e}")
    return uh, extra


def run_from_inp(inp_path: str, output_dir: str = "/app/output", debug: int = 0,
                 profile: str = None, units: str = None, parse_only: bool = False,
                 run_id: str = None, solver: str = "lu", rtol: float = 1e-10):
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
    run_id = run_id or str(uuid.uuid4())
    if profile_fc:
        write_run_status(output_dir, {
            "status": "started",
            "run_id": run_id,
            "solver_run": not parse_only,
            "solver": normalize_solver(solver),
            "rtol": normalize_rtol(rtol, normalize_solver(solver)),
            "note": "run em andamento; manifesto antigo neste diretorio nao e sucesso desta run",
        })

    model = parse_inp(inp_path, debug=debug, profile=profile)
    if debug:
        print(summary(model))
        print("---")

    if parse_only:
        payload = parse_report_payload(model, os.path.abspath(inp_path))
        payload["run_id"] = run_id
        payload["units"] = units_norm
        report_path = os.path.join(output_dir, "parse_report.json")
        with open(report_path, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, ensure_ascii=False)
            handle.write("\n")
        print(PARSE_NOT_SOLVED_BANNER)
        print("Parse report:", report_path)
        if profile_fc:
            write_run_status(output_dir, {
                "status": "parse_ok_not_solved",
                "run_id": run_id,
                "solver_run": False,
                "parse_report": report_path,
            })
        return {"parse_only": True, "model": model, "report": payload}

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
                model, step, output_base, debug, profile=profile,
                solver=solver, rtol=rtol,
            )
        results.append((result, extra))

    if profile_fc:
        _finalize_fc_run(model, output_dir, inp_path, units_norm, results, debug, run_id=run_id)

    return results


def _write_fc_linear_fields(domain, uh, v, sigma, epsilon, f_ext, a_form, model, output_base, debug,
                            constrained_rows=None, applied=None):
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
    # Simetria residual ja e aplicada em principals_algebraic.
    principals, vm = principals_algebraic(sigma_cells)
    prin_f = interpolate_cell_vector(domain, principals, "principal_stress")
    vm_f = interpolate_cell_scalar(domain, vm, "von_mises")

    cell_pvd = write_functions_vtk(
        domain.comm,
        [cauchy_f, strain_f, prin_f, vm_f],
        output_base + "_cell",
    )

    # R = K u - f_ext (assemble_vector esparso), somente nos DOFs prescritos (livres ~0).
    f_int = assemble_vector(fem.form(ufl.inner(sigma(uh), epsilon(v)) * ufl.dx))
    try:
        f_int.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
    except Exception:
        pass

    constrained_rows = constrained_rows or []
    applied = applied or {"Fx": 0.0, "Fy": 0.0, "Fz": 0.0}
    coords_by_nid = {nid: (n.x, n.y, n.z) for nid, n in model.nodes.items()}
    rx = residual_on_constrained(f_int, f_ext, constrained_rows, coords_by_nid)

    V = uh.function_space
    rf = fem.Function(V)
    rf.name = "reaction_force"
    rf.x.array[:] = 0.0
    mlen = min(rf.x.array.size, rx["rf_constrained"].size)
    rf.x.array[:mlen] = rx["rf_constrained"][:mlen]
    rf.x.scatter_forward()
    rf_pvd = write_functions_vtk(domain.comm, [rf], output_base + "_RF")

    totals = rx["totals"]
    sum_rf = {
        "Fx": totals["RFx"] + applied["Fx"],
        "Fy": totals["RFy"] + applied["Fy"],
        "Fz": totals["RFz"] + applied["Fz"],
    }
    # Momento das cargas aplicadas em torno da origem.
    mom_f = {"Mx": 0.0, "My": 0.0, "Mz": 0.0}
    for cl in (model.steps[0].loads if model.steps else []):
        if cl.node not in model.nodes:
            continue
        n = model.nodes[cl.node]
        fx = cl.value if cl.dof == 1 else 0.0
        fy = cl.value if cl.dof == 2 else 0.0
        fz = cl.value if cl.dof == 3 else 0.0
        mom_f["Mx"] += n.y * fz - n.z * fy
        mom_f["My"] += n.z * fx - n.x * fz
        mom_f["Mz"] += n.x * fy - n.y * fx
    mom_sum = {
        "Mx": rx["moment_R_origin"]["Mx"] + mom_f["Mx"],
        "My": rx["moment_R_origin"]["My"] + mom_f["My"],
        "Mz": rx["moment_R_origin"]["Mz"] + mom_f["Mz"],
    }
    force_char = 0.0
    moment_char = 0.0
    for cl in (model.steps[0].loads if model.steps else []):
        if cl.node not in model.nodes:
            continue
        n = model.nodes[cl.node]
        fx = cl.value if cl.dof == 1 else 0.0
        fy = cl.value if cl.dof == 2 else 0.0
        fz = cl.value if cl.dof == 3 else 0.0
        force_char += abs(fx) + abs(fy) + abs(fz)
        moment_char += abs(n.y * fz - n.z * fy) + abs(n.z * fx - n.x * fz) + abs(n.x * fy - n.y * fx)
    for ent in rx.get("entries") or []:
        val = float(ent.get("R", 0.0))
        comp = int(ent.get("component", 0))
        nid = int(ent.get("inp_node", 0))
        if nid not in model.nodes:
            continue
        n = model.nodes[nid]
        fx = val if comp == 0 else 0.0
        fy = val if comp == 1 else 0.0
        fz = val if comp == 2 else 0.0
        force_char += abs(val)
        moment_char += abs(n.y * fz - n.z * fy) + abs(n.z * fx - n.x * fz) + abs(n.x * fy - n.y * fx)

    react_path = output_base + "_reactions.json"
    payload = {
        "totals": totals,
        "applied_totals": applied,
        "sum_R_plus_F": sum_rf,
        "moment_R_origin": rx["moment_R_origin"],
        "moment_F_origin": mom_f,
        "moment_R_plus_F_origin": mom_sum,
        "n_constrained_dofs": rx["n_constrained_dofs"],
        "residual_free_dofs_max_abs": rx["residual_free_dofs_max_abs"],
        "force_char": force_char,
        "moment_char": moment_char,
        "entries": rx["entries"],
        "note": (
            "RF = (K_full u - f_ext) apenas em DOFs prescritos. "
            "Equilibrio de forca/momento: sum(R+F) e momentos em torno da origem ~0. "
            "Nao e aceite de roda/Luiz."
        ),
    }
    with open(react_path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")

    if debug:
        print("FC campos celula:", cell_pvd)
        print("FC RF:", totals, "R+F:", sum_rf, "free residual max:", rx["residual_free_dofs_max_abs"])

    honored = ["U", "S", "E", "PRINCIPAL", "VON_MISES", "RF"]
    artifacts = [
        {"path": output_base + ".pvd", "field": "U", "location": "nodal-CG1", "space": "Lagrange-P1"},
        {"path": cell_pvd, "field": "S,E,PRINCIPAL,VON_MISES", "location": "cell-DG0-constant-C3D4-P1",
         "tensor_order": "row-major XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ",
         "principals": "algebraic sigma1>=sigma2>=sigma3 with sign; tensor symmetrized before eig",
         "note": "nodal smoothed does not replace cell field"},
        {"path": rf_pvd, "field": "RF", "location": "nodal-constrained-dofs-only"},
        {"path": react_path, "field": "RF_totals", "location": "json"},
    ]
    return {
        "honored": honored,
        "artifacts": artifacts,
        "reactions": payload,
        "applied": applied,
        "stress_cells": sigma_cells,
        "principals": principals,
        "von_mises": vm,
    }


def _check_fc_force_balance(extra, rel=1e-5, abs_tol=1e-6):
    """Recusa sucesso se R+F (e momentos) nao fecharem no residual do proprio solver.

    Tolerancia numerica de consistencia (LU / residual), nao aceite Luiz/GateV.
    Forca/momento resultantes nulos (binario puro) usam a escala caracteristica
    das parcelas (|F_i| e |r x F_i| + |r x R_i|), nao o total que se cancela.
    """
    rx = extra.get("reactions") or {}
    sum_rf = rx.get("sum_R_plus_F") or {}
    mom = rx.get("moment_R_plus_F_origin") or {}
    applied = rx.get("applied_totals") or extra.get("applied") or {}
    totals = rx.get("totals") or {}
    scale_f = max(
        [
            abs_tol,
            abs(float(applied.get("Fx", 0.0))),
            abs(float(applied.get("Fy", 0.0))),
            abs(float(applied.get("Fz", 0.0))),
            abs(float(totals.get("RFx", 0.0))),
            abs(float(totals.get("RFy", 0.0))),
            abs(float(totals.get("RFz", 0.0))),
            abs(float(rx.get("force_char", 0.0))),
        ]
    )
    for key in ("Fx", "Fy", "Fz"):
        val = abs(float(sum_rf.get(key, 0.0)))
        if val > max(abs_tol, rel * scale_f):
            raise RuntimeError(
                "Equilibrio de forca nao fechou no perfil FC ({0}={1}, escala={2}). "
                "Nao e sucesso.".format(key, sum_rf.get(key), scale_f)
            )
    scale_m = max(
        [
            abs_tol,
            scale_f,  # fallback se nao houver binario
            abs(float(rx.get("moment_char", 0.0))),
            abs(float((rx.get("moment_F_origin") or {}).get("Mx", 0.0))),
            abs(float((rx.get("moment_F_origin") or {}).get("My", 0.0))),
            abs(float((rx.get("moment_F_origin") or {}).get("Mz", 0.0))),
            abs(float((rx.get("moment_R_origin") or {}).get("Mx", 0.0))),
            abs(float((rx.get("moment_R_origin") or {}).get("My", 0.0))),
            abs(float((rx.get("moment_R_origin") or {}).get("Mz", 0.0))),
        ]
    )
    for key in ("Mx", "My", "Mz"):
        val = abs(float(mom.get(key, 0.0)))
        if val > max(abs_tol, rel * scale_m):
            raise RuntimeError(
                "Equilibrio de momento nao fechou no perfil FC ({0}={1}, escala={2}). "
                "Nao e sucesso.".format(key, mom.get(key), scale_m)
            )


def _stack_versions():
    info = {"python": sys.version.split()[0]}
    for name in ("dolfinx", "ufl", "basix", "ffcx", "petsc4py", "mpi4py", "numpy"):
        try:
            mod = __import__(name)
            info[name] = getattr(mod, "__version__", "?")
        except Exception:
            info[name] = "ausente"
    return info


def _finalize_fc_run(model, output_dir, inp_path, units_norm, results, debug, run_id=None):
    rec = model.fc_record
    honored = []
    artifacts = []
    reactions = None
    solver_info = None
    for _res, extra in results:
        honored = list(extra.get("honored") or [])
        artifacts.extend(extra.get("artifacts") or [])
        reactions = extra.get("reactions") or reactions
        solver_info = extra.get("solver_info") or solver_info

    requested = required_output_fields(rec) if rec else ["U", "S", "PRINCIPAL", "VON_MISES"]
    check_honored_outputs(requested, honored)

    payload = {
        "status": "ok",
        "run_id": run_id,
        "solver_run": True,
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
            "S": "Cauchy cell DG0, constant on C3D4 P1; not nodal-smoothed; symmetrized before principals",
            "E": "infinitesimal strain cell DG0",
            "PRINCIPAL": "algebraic sigma1>=sigma2>=sigma3 with sign, cell DG0",
            "VON_MISES": "cell DG0",
            "RF": "residual K_full u - f_ext at constrained DOFs only, plus force/moment totals",
        },
        "static_increments": rec.static_increments if rec else None,
        "reactions": reactions,
        "solver": solver_info,
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
            "independent_benchmark": "axial+compressao analitica em fixture sintetica; GateV pendente",
            "material_plan_70000": "nao corrigido automaticamente se o .inp trouxer outro E",
            "mpi": "serial only (MPI size=1)",
        },
    }
    manifest_path = write_run_manifest(output_dir, payload)
    required = [os.path.basename(a["path"]) for a in artifacts if a.get("path")]
    verify_run_artifacts(manifest_path, required)
    status_body = {
        "status": "ok",
        "run_id": run_id,
        "solver_run": True,
        "manifest": manifest_path,
    }
    if solver_info:
        status_body["solver"] = solver_info
    write_run_status(output_dir, status_body)
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
    output_dir = resolve_output_dir(args.output_dir)
    run_id = str(uuid.uuid4())
    profile_fc = (args.profile or "").strip().lower() == "fc"
    try:
        run_from_inp(
            args.inp_file,
            output_dir=output_dir,
            debug=args.debug,
            profile=args.profile,
            units=args.units,
            parse_only=args.parse_only,
            run_id=run_id,
            solver=args.solver,
            rtol=args.rtol,
        )
    except FCProfileError as exc:
        print("ERRO FC: {0}".format(exc))
        if profile_fc:
            try:
                write_run_status(prepare_output_dir(output_dir), {
                    "status": "failed",
                    "run_id": run_id,
                    "solver_run": False if args.parse_only else True,
                    "error": str(exc),
                    "error_class": "FCProfileError",
                    "note": "falha: manifesto de sucesso anterior neste diretorio nao vale para esta run",
                })
            except Exception:
                pass
        sys.exit(2)
    except RuntimeError as exc:
        print("ERRO: {0}".format(exc))
        if profile_fc:
            try:
                failed = {
                    "status": "failed",
                    "run_id": run_id,
                    "solver_run": True,
                    "error": str(exc),
                    "error_class": "RuntimeError",
                    "note": "falha: manifesto de sucesso anterior neste diretorio nao vale para esta run",
                }
                if getattr(exc, "solver_info", None):
                    failed["solver"] = exc.solver_info
                else:
                    failed["solver"] = {"type": getattr(args, "solver", None), "rtol": getattr(args, "rtol", None)}
                write_run_status(prepare_output_dir(output_dir), failed)
            except Exception:
                pass
        sys.exit(3)
    except Exception as exc:
        print("ERRO: {0}: {1}".format(type(exc).__name__, exc))
        if profile_fc:
            try:
                write_run_status(prepare_output_dir(output_dir), {
                    "status": "failed",
                    "run_id": run_id,
                    "solver_run": True,
                    "error": "{0}: {1}".format(type(exc).__name__, exc),
                    "error_class": type(exc).__name__,
                    "note": "falha: manifesto de sucesso anterior neste diretorio nao vale para esta run",
                })
            except Exception:
                pass
        sys.exit(3)
