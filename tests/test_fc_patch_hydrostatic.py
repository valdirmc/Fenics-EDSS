"""T0 hidrostatico e T2 patch test do solver FC (C3D4 linear).

Os testes passam pelo leitor .inp do perfil FC (o motor nao gera INP de
producao; fixtures de verificacao em tests/fixtures/ sao permitidas).

Tolerancias dos asserts: as formulas C4/C5/C6 propostas por Luiz (13/11),
com u_ref/sigma_ref PROPOSTOS abaixo, ainda pendentes de congelamento pelo
Valdir (D8). Nao sao mais folgadas que essas formulas. NAO e aceite oficial,
NAO congela criterio, NAO aprova Gate V.

COORD (1e-9*L_ref + 2*ulp) e proposta pendente, documentada, nao assertada
como criterio congelado.
"""

from __future__ import annotations

import json
import os
import sys

import pytest

pytest.importorskip("dolfinx")
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fc_fields import principals_algebraic  # noqa: E402
from fc_verify_meshes import (  # noqa: E402
    E,
    L,
    NU,
    PRESSURE,
    even_permutations_4,
    hydrostatic_alpha,
    hydrostatic_eps_vol,
    hydrostatic_u,
    make_mesh,
    patch_sigma,
    patch_strain,
    patch_u,
    permute_cells,
    tet_volumes,
    build_hydrostatic,
    build_patch_dirichlet,
    build_patch_traction,
)
from run_from_inp import run_from_inp  # noqa: E402

# Equilibrio de forca: consistencia numerica, nao e C0-C11.
ABS_F = 1e-6
FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")

# PROPOSTA Luiz 13/11 (revisao da sprint 2), pendente de decisao do Valdir (D8).
# T2: u_ref = max_k |u_ex,k| = 0,222 mm; sigma_ref = max_ij |sigma_ex,ij| = 289,473684 MPa (syz).
# T0: u_ref = p L (1-2 nu)/E = 4,857e-4 mm; sigma_ref = p = 10 MPa.
T2_U_REF = 0.222  # mm, PROPOSTA
T2_SIGMA_REF = 289.473684  # MPa, PROPOSTA
T0_U_REF = 4.857e-4  # mm, PROPOSTA
T0_SIGMA_REF = 10.0  # MPa, PROPOSTA
T0_P = PRESSURE


def c4_tol(u_ref, u_ref_comp):
    """C4 proposto: |du_k| <= 1e-6*u_ref + 1e-5*|u_ref,k|. Nao congelado."""
    return 1e-6 * float(u_ref) + 1e-5 * np.abs(np.asarray(u_ref_comp, dtype=float))


def c5_tol(sigma_ref, s_ref_comp):
    """C5 proposto: |dS_ij| <= 1e-5*sigma_ref + 1e-4*|S_ij,ref|. Nao congelado."""
    return 1e-5 * float(sigma_ref) + 1e-4 * np.abs(np.asarray(s_ref_comp, dtype=float))


def c6_tol_invariants(sigma_ref, principals_ref):
    """C6 proposto: |dS_k| <= 1e-5*sigma_ref + 1e-4*max(|S1|,|S3|)_ref. Nao congelado."""
    pr = np.asarray(principals_ref, dtype=float)
    scale = np.maximum(np.abs(pr[..., 0]), np.abs(pr[..., 2]))
    return 1e-5 * float(sigma_ref) + 1e-4 * scale


def c6_vm_hydro_tol(p):
    """C6/T0 proposto no pos-processador e na solucao em memoria: sigma_vm <= 1e-12*|p|."""
    return 1e-12 * abs(float(p))


def assert_c4(du, u_ref, u_b, label):
    du = np.asarray(du, dtype=float)
    tol = c4_tol(u_ref, u_b)
    worst = float(np.max(np.abs(du) - tol))
    assert np.all(np.abs(du) <= tol), "{0}: C4 proposto violado (max(|du|-tol)={1:.6e})".format(label, worst)


def assert_c5(dS, sigma_ref, s_b, label):
    dS = np.asarray(dS, dtype=float)
    tol = c5_tol(sigma_ref, s_b)
    worst = float(np.max(np.abs(dS) - tol))
    assert np.all(np.abs(dS) <= tol), "{0}: C5 proposto violado (max(|dS|-tol)={1:.6e})".format(label, worst)


def assert_c6_invariants(dpr, dvm, sigma_ref, pr_b, label):
    dpr = np.asarray(dpr, dtype=float)
    pr_b = np.asarray(pr_b, dtype=float)
    tol = np.asarray(c6_tol_invariants(sigma_ref, pr_b), dtype=float)
    if dpr.ndim == 2 and tol.ndim == 1:
        tol_cmp = tol[:, None]
    else:
        tol_cmp = tol
    worst = float(np.max(np.abs(dpr) - tol_cmp))
    assert np.all(np.abs(dpr) <= tol_cmp), (
        "{0}: C6 proposto (principais) violado (max(|dpr|-tol)={1:.6e})".format(label, worst)
    )
    vm_tol = float(np.max(tol))
    assert float(np.max(np.abs(np.asarray(dvm, dtype=float)))) <= vm_tol, (
        "{0}: C6 proposto (VM) violado".format(label)
    )


def read_fixture(name: str) -> str:
    path = os.path.join(FIXTURES, name)
    with open(path, encoding="utf-8") as handle:
        return handle.read()


def _run(tmp_path, text, name, solver="lu", rtol=1e-10):
    inp = tmp_path / name
    inp.write_text(text, encoding="utf-8")
    out = tmp_path / ("out_" + name.replace(".inp", ""))
    out.mkdir()
    results = run_from_inp(
        str(inp),
        output_dir=str(out),
        debug=0,
        profile="fc",
        units="mm-N-MPa-t-s",
        solver=solver,
        rtol=rtol,
    )
    return results, out


def _u_and_coords(uh):
    coords = np.asarray(uh.function_space.mesh.geometry.x, dtype=float)
    uarr = np.asarray(uh.x.array, dtype=float).reshape((-1, 3))
    if coords.shape[0] != uarr.shape[0]:
        raise RuntimeError("Nos da geometria ({0}) != dofs de U ({1})".format(coords.shape[0], uarr.shape[0]))
    return coords, uarr


def _match_original(coords_mesh, coords_orig, atol=1e-9):
    """Permutacao da malha DOLFINx de volta aos nos originais do INP."""
    index = np.empty(coords_orig.shape[0], dtype=int)
    used = np.zeros(coords_mesh.shape[0], dtype=bool)
    for i, x in enumerate(coords_orig):
        dist = np.linalg.norm(coords_mesh - x, axis=1)
        dist[used] = np.inf
        j = int(np.argmin(dist))
        if not np.isfinite(dist[j]) or dist[j] > atol:
            raise RuntimeError("No original {0} nao encontrado (min dist={1})".format(i, dist[j]))
        used[j] = True
        index[i] = j
    if not used.all():
        raise RuntimeError("Alguns vertices da malha nao casaram com o INP.")
    return index


def _strain_vol_from_u(coords, cells, u_orig):
    ev = np.empty(cells.shape[0], dtype=float)
    for i, cell in enumerate(cells):
        x = coords[cell]
        u = u_orig[cell]
        jac = np.column_stack((x[1] - x[0], x[2] - x[0], x[3] - x[0]))
        grel = np.column_stack((u[1] - u[0], u[2] - u[0], u[3] - u[0]))
        grad = grel @ np.linalg.inv(jac)
        ev[i] = float(np.trace(0.5 * (grad + grad.T)))
    return ev


def _metrics_common(extra, out):
    rx = extra.get("reactions") or {}
    info = extra.get("solver_info") or {}
    sum_rf = rx.get("sum_R_plus_F") or {}
    mom = rx.get("moment_R_plus_F_origin") or {}
    principals = np.asarray(extra["principals"], dtype=float)
    vm = np.asarray(extra["von_mises"], dtype=float)
    sigma = np.asarray(extra["stress_cells"], dtype=float)
    status = json.loads((out / "run_status.json").read_text(encoding="utf-8"))
    manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    finite = bool(np.isfinite(sigma).all() and np.isfinite(principals).all() and np.isfinite(vm).all())
    return {
        "sum_R_plus_F": {k: float(sum_rf.get(k, 0.0)) for k in ("Fx", "Fy", "Fz")},
        "sum_R_plus_F_max_abs": float(max(abs(float(sum_rf.get(k, 0.0))) for k in ("Fx", "Fy", "Fz"))),
        "moment_R_plus_F_origin": {k: float(mom.get(k, 0.0)) for k in ("Mx", "My", "Mz")},
        "moment_max_abs": float(max(abs(float(mom.get(k, 0.0))) for k in ("Mx", "My", "Mz"))),
        "residual_free_dofs_max_abs": float(rx.get("residual_free_dofs_max_abs") or 0.0),
        "true_residual_rel": info.get("true_residual_rel"),
        "true_residual_meets_rtol": info.get("true_residual_meets_rtol"),
        "solver_type": info.get("type"),
        "solver_iterations": info.get("iterations"),
        "n_cells": int(sigma.shape[0]),
        "fields_finite": finite,
        "status": status.get("status"),
        "run_id": status.get("run_id"),
        "units": manifest.get("units"),
        "principals_ordered": bool(np.all(principals[:, 0] >= principals[:, 1] - 1e-14)
                                   and np.all(principals[:, 1] >= principals[:, 2] - 1e-14)),
    }


def _assert_c11_and_balance(m, abs_f=ABS_F):
    assert m["status"] == "ok"
    assert m["units"] == "mm-N-MPa-t-s"
    assert m["fields_finite"] is True
    assert m["n_cells"] > 0
    assert m["principals_ordered"] is True
    assert m["sum_R_plus_F_max_abs"] == pytest.approx(0.0, abs=abs_f)
    # Momento: o hidrostatico cancela parcelas grandes; folga um pouco maior.
    assert m["moment_max_abs"] == pytest.approx(0.0, abs=1e-5)


def test_t0_posprocessador_hidrostatico_proposto():
    """Ja coberto em test_fc_profile; reforco do criterio proposto VM <= 1e-12*|p|."""
    for pval in (PRESSURE, -PRESSURE, 1.0, 100.0):
        sigma = np.eye(3, dtype=float)[None, :, :] * (-pval if pval == PRESSURE else pval)
        # Para p de pressao positiva, sigma = -p I. Mantem os dois sinais.
        if pval in (PRESSURE,):
            sigma = -pval * np.eye(3, dtype=float)[None, :, :]
        else:
            sigma = pval * np.eye(3, dtype=float)[None, :, :]
        pr, vm = principals_algebraic(sigma)
        scale = abs(float(pval))
        assert pr[0, 0] == pytest.approx(sigma[0, 0, 0])
        assert pr[0, 2] == pytest.approx(sigma[0, 0, 0])
        assert float(vm[0]) <= 1e-12 * scale


def test_t2_patch_dirichlet_contorno(tmp_path):
    pack = build_patch_dirichlet()
    results, out = _run(tmp_path, pack["text"], "t2_patch_block.inp")
    uh, extra = results[0]
    coords_mesh, u_mesh = _u_and_coords(uh)
    index = _match_original(coords_mesh, pack["coords"])
    u_orig = u_mesh[index]
    u_ex = patch_u(pack["coords"])
    du = u_orig - u_ex
    max_du = float(np.max(np.abs(du)))
    max_du_interior = float(np.max(np.abs(du[~pack["boundary"]]))) if np.any(~pack["boundary"]) else 0.0
    sigma = np.asarray(extra["stress_cells"], dtype=float)
    # Tensao DG0 pode seguir a ordem das celulas da malha, nao do INP.
    # Como o estado e constante, basta comparar todas as celulas ao analitico.
    s_ex = patch_sigma()
    dS = sigma - s_ex[None, :, :]
    max_dS = float(np.max(np.abs(dS)))
    max_dS_named = {
        "sxx": float(np.max(np.abs(dS[:, 0, 0]))),
        "syy": float(np.max(np.abs(dS[:, 1, 1]))),
        "szz": float(np.max(np.abs(dS[:, 2, 2]))),
        "sxy": float(np.max(np.abs(dS[:, 0, 1]))),
        "sxz": float(np.max(np.abs(dS[:, 0, 2]))),
        "syz": float(np.max(np.abs(dS[:, 1, 2]))),
    }
    pr = np.asarray(extra["principals"], dtype=float)
    pr_ex, vm_ex = principals_algebraic(s_ex[None, :, :])
    max_dpr = float(np.max(np.abs(pr - pr_ex[0][None, :])))
    max_dvm = float(np.max(np.abs(np.asarray(extra["von_mises"]) - vm_ex[0])))
    m = _metrics_common(extra, out)
    _assert_c11_and_balance(m)
    assert pack["meta"]["n_elements"] == 384
    assert int(np.count_nonzero(~pack["boundary"])) >= 27
    assert pack["vols"].min() > 0.0
    u_ex = patch_u(pack["coords"])
    assert_c4(du, T2_U_REF, u_ex, "T2 Dirichlet")
    assert_c4(du[~pack["boundary"]], T2_U_REF, u_ex[~pack["boundary"]], "T2 Dirichlet interiores")
    assert_c5(dS, T2_SIGMA_REF, np.broadcast_to(s_ex, dS.shape), "T2 Dirichlet")
    pr_ex, vm_ex = principals_algebraic(s_ex[None, :, :])
    assert_c6_invariants(pr - pr_ex[0][None, :], np.asarray(extra["von_mises"]) - vm_ex[0],
                         T2_SIGMA_REF, np.broadcast_to(pr_ex[0], pr.shape), "T2 Dirichlet")
    # Cisalhamentos distintos: se houvesse troca de componentes, dS explodiria.
    shear = patch_strain()
    assert len({round(shear[0, 1], 12), round(shear[0, 2], 12), round(shear[1, 2], 12)}) == 3
    assert abs(shear[0, 1]) > 0 and abs(shear[0, 2]) > 0 and abs(shear[1, 2]) > 0
    print(
        "T2 PATCH DIRICHLET MEASURED: n_nodes={0} n_elem={1} n_interior={2} "
        "min_vol={3:.6e} max|du|={4:.6e} max|du|_int={5:.6e} max|dS|={6:.6e} "
        "dS_named={7} max|dS1-3|={8:.6e} max|dVM|={9:.6e} "
        "sum|R+F|={10:.6e} mom={11:.6e} res_free={12:.6e}".format(
            pack["meta"]["n_nodes"], pack["meta"]["n_elements"],
            int(np.count_nonzero(~pack["boundary"])), pack["vols"].min(),
            max_du, max_du_interior, max_dS, max_dS_named, max_dpr, max_dvm,
            m["sum_R_plus_F_max_abs"], m["moment_max_abs"],
            m["residual_free_dofs_max_abs"],
        )
    )


def test_t2_patch_12_permutacoes_pares(tmp_path):
    pack = build_patch_dirichlet()
    s_ex = patch_sigma()
    reports = []
    for k, perm in enumerate(even_permutations_4()):
        cells = permute_cells(pack["cells"], perm)
        vols = tet_volumes(pack["coords"], cells)
        assert np.all(vols > 0.0), "perm {0} volume <= 0: {1}".format(perm, vols.min())
        # Reconstroi o texto com a conectividade permutada, mesmos BC.
        from fc_verify_meshes import inp_text

        text = inp_text(
            pack["coords"], cells, pack["bcs"], {}, pack["meta"],
            heading="T2 patch C3D4, permutacao par {0} {1}".format(k, perm),
        )
        results, out = _run(tmp_path, text, "t2_perm_{0}.inp".format(k))
        uh, extra = results[0]
        coords_mesh, u_mesh = _u_and_coords(uh)
        index = _match_original(coords_mesh, pack["coords"])
        du = u_mesh[index] - patch_u(pack["coords"])
        dS = np.asarray(extra["stress_cells"]) - s_ex[None, :, :]
        rec = {
            "perm": list(perm),
            "max_du": float(np.max(np.abs(du))),
            "max_dS": float(np.max(np.abs(dS))),
            "min_vol": float(vols.min()),
        }
        reports.append(rec)
        assert_c4(du, T2_U_REF, patch_u(pack["coords"]), "T2 perm {0}".format(perm))
        assert_c5(dS, T2_SIGMA_REF, np.broadcast_to(s_ex, dS.shape), "T2 perm {0}".format(perm))
    worst_u = max(r["max_du"] for r in reports)
    worst_s = max(r["max_dS"] for r in reports)
    print(
        "T2 PATCH 12 PERMS MEASURED: n_perm=12 worst|du|={0:.6e} worst|dS|={1:.6e} reports={2}".format(
            worst_u, worst_s, reports
        )
    )


def test_t2_patch_tracao_dirichlet_parcial(tmp_path):
    pack = build_patch_traction()
    text = read_fixture("t2_patch_tracao.inp")
    results, out = _run(tmp_path, text, "t2_patch_tracao.inp")
    uh, extra = results[0]
    coords_mesh, u_mesh = _u_and_coords(uh)
    index = _match_original(coords_mesh, pack["coords"])
    du = u_mesh[index] - patch_u(pack["coords"])
    free = ~pack["dirichlet"]
    max_du = float(np.max(np.abs(du)))
    max_du_free = float(np.max(np.abs(du[free]))) if np.any(free) else 0.0
    dS = np.asarray(extra["stress_cells"]) - patch_sigma()[None, :, :]
    max_dS = float(np.max(np.abs(dS)))
    m = _metrics_common(extra, out)
    _assert_c11_and_balance(m)
    u_ex = patch_u(pack["coords"])
    assert_c4(du, T2_U_REF, u_ex, "T2 tracao")
    assert_c4(du[free], T2_U_REF, u_ex[free], "T2 tracao livres")
    assert_c5(dS, T2_SIGMA_REF, np.broadcast_to(patch_sigma(), dS.shape), "T2 tracao")
    print(
        "T2 PATCH TRACAO MEASURED: n_elem={0} n_free={1} max|du|={2:.6e} "
        "max|du|_free={3:.6e} max|dS|={4:.6e} sum|R+F|={5:.6e} mom={6:.6e} res_free={7:.6e}".format(
            pack["meta"]["n_elements"], int(np.count_nonzero(free)),
            max_du, max_du_free, max_dS,
            m["sum_R_plus_F_max_abs"], m["moment_max_abs"],
            m["residual_free_dofs_max_abs"],
        )
    )


def _hydro_check(tmp_path, solver="lu", name="t0_hydrostatic_block.inp"):
    pack = build_hydrostatic()
    results, out = _run(tmp_path, pack["text"], name, solver=solver)
    uh, extra = results[0]
    p = pack["pressure"]
    s_ex = -p * np.eye(3, dtype=float)
    sigma = np.asarray(extra["stress_cells"], dtype=float)
    dS = sigma - s_ex[None, :, :]
    max_dS = float(np.max(np.abs(dS)))
    vm = np.asarray(extra["von_mises"], dtype=float)
    max_vm = float(np.max(np.abs(vm)))
    pr = np.asarray(extra["principals"], dtype=float)
    max_dpr = float(np.max(np.abs(pr + p)))
    coords_mesh, u_mesh = _u_and_coords(uh)
    index = _match_original(coords_mesh, pack["coords"])
    du = u_mesh[index] - hydrostatic_u(pack["coords"], p)
    max_du = float(np.max(np.abs(du)))
    ev = _strain_vol_from_u(pack["coords"], pack["cells"], u_mesh[index])
    ev_ex = hydrostatic_eps_vol(p)
    max_dev = float(np.max(np.abs(ev - ev_ex)))
    # Tambem via lei de Hooke inversa: tr(eps) = (1-2 nu)/E * tr(sigma)
    ev_from_s = ((1.0 - 2.0 * NU) / E) * np.trace(sigma, axis1=1, axis2=2)
    max_dev_s = float(np.max(np.abs(ev_from_s - ev_ex)))
    m = _metrics_common(extra, out)
    _assert_c11_and_balance(m)
    assert pack["vols"].min() > 0.0
    u_ex = hydrostatic_u(pack["coords"], p)
    assert_c4(du, T0_U_REF, u_ex, "T0 {0}".format(solver))
    assert_c5(dS, T0_SIGMA_REF, np.broadcast_to(s_ex, dS.shape), "T0 {0}".format(solver))
    pr_ex, vm_ex = principals_algebraic(s_ex[None, :, :])
    assert_c6_invariants(pr - pr_ex[0][None, :], vm - vm_ex[0],
                         T0_SIGMA_REF, np.broadcast_to(pr_ex[0], pr.shape), "T0 {0}".format(solver))
    luiz_vm_tol = c6_vm_hydro_tol(p)
    # C6/T0 proposto: sigma_vm em memoria <= 1e-12*|p| (nao julga o CCX impresso).
    assert max_vm <= luiz_vm_tol, "T0 {0}: C6 proposto VM {1:.6e} > {2:.6e}".format(solver, max_vm, luiz_vm_tol)
    meets_luiz_vm = bool(max_vm <= luiz_vm_tol)
    # Deformacao volumetrica: consistencia do campo no espaco P1 (nao e C0-C11).
    assert max_dev <= 1e-12
    report = {
        "solver": solver,
        "p": p,
        "n_nodes": pack["meta"]["n_nodes"],
        "n_elements": pack["meta"]["n_elements"],
        "min_volume": pack["vols"].min(),
        "alpha": hydrostatic_alpha(p),
        "eps_vol_exact": ev_ex,
        "max_du": max_du,
        "max_dS": max_dS,
        "max_vm": max_vm,
        "max_d_principals": max_dpr,
        "max_d_eps_vol_u": max_dev,
        "max_d_eps_vol_sigma": max_dev_s,
        "luiz_T0_vm_tol": luiz_vm_tol,
        "meets_luiz_T0_vm_1e12_p": meets_luiz_vm,
        **m,
    }
    print("T0 HIDROSTATICO MEASURED ({0}): {1}".format(solver, json.dumps(report, sort_keys=True)))
    return report


def test_t0_hidrostatico_solver_lu(tmp_path):
    _hydro_check(tmp_path, solver="lu", name="t0_hydro_lu.inp")


def test_t0_hidrostatico_solver_cg_gamg(tmp_path):
    report = _hydro_check(tmp_path, solver="cg-gamg", name="t0_hydro_cg.inp")
    assert report["solver_type"] == "cg-gamg"
    assert report["true_residual_rel"] is not None
    assert report["true_residual_rel"] <= 1e-10
    assert report["true_residual_meets_rtol"] is True


def test_fixtures_inp_existem_e_sao_c3d4():
    needed = (
        "t2_patch_block.inp",
        "t2_patch_tracao.inp",
        "t0_hydrostatic_block.inp",
    )
    builders = {
        "t2_patch_block.inp": build_patch_dirichlet,
        "t2_patch_tracao.inp": build_patch_traction,
        "t0_hydrostatic_block.inp": build_hydrostatic,
    }
    for name in needed:
        path = os.path.join(FIXTURES, name)
        assert os.path.isfile(path), path
        text = read_fixture(name)
        assert "TYPE=C3D4" in text
        assert "*DLOAD" not in text
        assert "C3D8" not in text and "C3D10" not in text
        generated = builders[name]()["text"]
        assert text == generated, "{0} diverge do gerador atual".format(name)


def test_t2_e_t0_via_fixture_arquivo(tmp_path):
    """Resolve e carrega as fixtures gravadas, inclusive a de tracao."""
    cases = (
        ("t2_patch_block.inp", "fix_patch.inp", "patch"),
        ("t2_patch_tracao.inp", "fix_tracao.inp", "tracao"),
        ("t0_hydrostatic_block.inp", "fix_hydro.inp", "hydro"),
    )
    for src_name, name, kind in cases:
        src = os.path.join(FIXTURES, src_name)
        assert os.path.isfile(src), src
        text = read_fixture(src_name)
        results, out = _run(tmp_path, text, name)
        _uh, extra = results[0]
        m = _metrics_common(extra, out)
        _assert_c11_and_balance(m)
        if kind in {"patch", "tracao"}:
            dS = np.asarray(extra["stress_cells"]) - patch_sigma()[None, :, :]
            assert_c5(dS, T2_SIGMA_REF, np.broadcast_to(patch_sigma(), dS.shape), kind)
        else:
            vm = np.asarray(extra["von_mises"])
            assert float(np.max(np.abs(vm))) <= c6_vm_hydro_tol(T0_P)
