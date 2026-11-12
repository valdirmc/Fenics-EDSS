"""Solver FC: prisma C3D4 com forcas externas e DOFs livres.

Referencia: lei de Hooke isotropica, estado uniaxial.
Cargas nodais consistentes em triangulo linear sob tracao constante:
F_i = A * t / 3 (Zienkiewicz & Taylor, The Finite Element Method,
vetor de carga consistente do triangulo linear / CST).

Campo de deslocamento linear (e portanto exato em P1/C3D4):
  sigma_zz = F/A, demais componentes nulas
  eps_zz = sigma_zz / E
  eps_xx = eps_yy = -nu * sigma_zz / E
  u_z = eps_zz * z, u_x = eps_xx * x, u_y = eps_yy * y

Tolerancia numerica dos testes (LU direto, campo no espaco FE):
  abs=1e-6 (MPa ou mm) ou rel=1e-8. NAO e aceite de Luiz nem GateV aprovado.

Malha sintetica em string; nao toca o INP da roda nem exporta INP de produto.
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

from inp_parser import parse_inp  # noqa: E402
from run_from_inp import run_from_inp  # noqa: E402


E = 70000.0
NU = 0.33
AREA = 0.5
LENGTH = 1.0
FZ = -50.0
SIGMA_ZZ = FZ / AREA  # -100 MPa
EPS_ZZ = SIGMA_ZZ / E
EPS_LAT = -NU * SIGMA_ZZ / E
SHARE = FZ / 3.0

# Tolerancia justificada: solucao no espaco P1, fatoracao LU.
ABS_U = 1e-8
ABS_S = 1e-6
ABS_F = 1e-6


def _prism_inp(*, with_cload=True, split_cload=False, prescribed_top=False, with_bc=True, with_output=True):
    lines = [
        "*NODE, NSET=Nall",
        "1, 0, 0, 0",
        "2, 1, 0, 0",
        "3, 0, 1, 0",
        "4, 0, 0, 1",
        "5, 1, 0, 1",
        "6, 0, 1, 1",
        "*ELEMENT, TYPE=C3D4, ELSET=Eall",
        "1, 1, 2, 3, 4",
        "2, 2, 3, 4, 5",
        "3, 3, 4, 5, 6",
        "*NSET, NSET=Bottom",
        "1, 2, 3",
        "*NSET, NSET=Top",
        "4, 5, 6",
        "*MATERIAL, NAME=Mat",
        "*ELASTIC",
        "{0}, {1}".format(E, NU),
        "*SOLID SECTION, ELSET=Eall, MATERIAL=Mat",
        "*STEP, INC=100",
        "*STATIC",
        "1.0, 1.0",
    ]
    if with_bc:
        lines.extend([
            "*BOUNDARY",
            "1, 1, 1, 0.0",
            "1, 2, 2, 0.0",
            "1, 3, 3, 0.0",
            "2, 2, 2, 0.0",
            "2, 3, 3, 0.0",
            "3, 3, 3, 0.0",
        ])
        if prescribed_top:
            uz = EPS_ZZ * LENGTH
            lines.extend([
                "4, 3, 3, {0}".format(uz),
                "5, 3, 3, {0}".format(uz),
                "6, 3, 3, {0}".format(uz),
            ])
    if with_cload:
        lines.append("*CLOAD")
        if split_cload:
            half = SHARE / 2.0
            lines.append("4, 3, {0}".format(half))
            lines.append("4, 3, {0}".format(SHARE - half))
        else:
            lines.append("4, 3, {0}".format(SHARE))
        lines.append("5, 3, {0}".format(SHARE))
        lines.append("6, 3, {0}".format(SHARE))
    if with_output:
        lines.extend([
            "*OUTPUT, FREQUENCY=1",
            "*NODE FILE",
            "U, RF",
            "*EL FILE",
            "S, E",
            "*NODE PRINT, NSET=Bottom, TOTALS=ONLY",
            "RF",
        ])
    lines.append("*END STEP")
    return "\n".join(lines) + "\n"


def _run(tmp_path, text, name="prism.inp", parse_only=False, solver="lu", rtol=1e-10):
    inp = tmp_path / name
    inp.write_text(text, encoding="utf-8")
    out = tmp_path / ("out_" + name)
    out.mkdir()
    return run_from_inp(
        str(inp),
        output_dir=str(out),
        debug=1,
        profile="fc",
        units="mm-N-MPa-t-s",
        parse_only=parse_only,
        solver=solver,
        rtol=rtol,
    ), out


def test_axial_compressao_cload_dofs_livres(tmp_path):
    results, out = _run(tmp_path, _prism_inp(split_cload=True), "axial.inp")
    uh, extra = results[0]
    # Deslocamentos nos nos INP via valores da funcao (ordem do espaco).
    coords = uh.function_space.mesh.geometry.x
    uarr = uh.x.array.reshape((-1, 3))
    # Localiza nos pelos coords da malha
    def u_at(xyz):
        dist = np.linalg.norm(coords - np.array(xyz), axis=1)
        i = int(np.argmin(dist))
        assert dist[i] < 1e-9
        return uarr[i]

    # Top node 4 (0,0,1): u = (0, 0, eps_zz)
    u4 = u_at((0.0, 0.0, 1.0))
    u5 = u_at((1.0, 0.0, 1.0))
    u6 = u_at((0.0, 1.0, 1.0))
    assert u4[2] == pytest.approx(EPS_ZZ * LENGTH, abs=ABS_U)
    assert u5[0] == pytest.approx(EPS_LAT * 1.0, abs=ABS_U)
    assert u5[2] == pytest.approx(EPS_ZZ * LENGTH, abs=ABS_U)
    assert u6[1] == pytest.approx(EPS_LAT * 1.0, abs=ABS_U)
    assert u6[2] == pytest.approx(EPS_ZZ * LENGTH, abs=ABS_U)
    u2 = u_at((1.0, 0.0, 0.0))
    assert u2[0] == pytest.approx(EPS_LAT, abs=ABS_U)
    assert u2[2] == pytest.approx(0.0, abs=ABS_U)

    sigma = extra["stress_cells"]
    principals = extra["principals"]
    vm = extra["von_mises"]
    assert sigma.shape[0] == 3
    for c in range(3):
        S = sigma[c]
        assert S[2, 2] == pytest.approx(SIGMA_ZZ, abs=ABS_S)
        assert S[0, 0] == pytest.approx(0.0, abs=ABS_S)
        assert S[1, 1] == pytest.approx(0.0, abs=ABS_S)
        assert S[0, 1] == pytest.approx(0.0, abs=ABS_S)
        assert S[0, 2] == pytest.approx(0.0, abs=ABS_S)
        assert S[1, 2] == pytest.approx(0.0, abs=ABS_S)
        # simetria
        assert abs(S[0, 1] - S[1, 0]) < ABS_S
        assert principals[c, 0] == pytest.approx(0.0, abs=ABS_S)
        assert principals[c, 1] == pytest.approx(0.0, abs=ABS_S)
        assert principals[c, 2] == pytest.approx(SIGMA_ZZ, abs=ABS_S)
        assert principals[c, 0] >= principals[c, 1] >= principals[c, 2]
        assert vm[c] == pytest.approx(abs(SIGMA_ZZ), abs=ABS_S)

    rx = extra["reactions"]
    assert rx["applied_totals"]["Fz"] == pytest.approx(FZ, abs=ABS_F)
    assert rx["totals"]["RFz"] == pytest.approx(-FZ, abs=ABS_F)
    assert rx["totals"]["RFx"] == pytest.approx(0.0, abs=ABS_F)
    assert rx["totals"]["RFy"] == pytest.approx(0.0, abs=ABS_F)
    assert rx["sum_R_plus_F"]["Fz"] == pytest.approx(0.0, abs=ABS_F)
    assert rx["sum_R_plus_F"]["Fx"] == pytest.approx(0.0, abs=ABS_F)
    assert rx["sum_R_plus_F"]["Fy"] == pytest.approx(0.0, abs=ABS_F)
    assert rx["moment_R_plus_F_origin"]["Mx"] == pytest.approx(0.0, abs=1e-5)
    assert rx["moment_R_plus_F_origin"]["My"] == pytest.approx(0.0, abs=1e-5)
    assert rx["moment_R_plus_F_origin"]["Mz"] == pytest.approx(0.0, abs=1e-5)
    # DOFs livres: residual deve ser numerico pequeno (nao o equilibrio trivial todos-fixos)
    assert rx["residual_free_dofs_max_abs"] < 1e-5

    manifest = json.loads((out / "run_manifest.json").read_text())
    status = json.loads((out / "run_status.json").read_text())
    assert status["status"] == "ok"
    assert status["run_id"] == manifest["run_id"]
    assert manifest["units"] == "mm-N-MPa-t-s"
    assert manifest["materials"]["Mat"]["E"] == E


def test_axial_deslocamento_prescrito_nao_zero(tmp_path):
    results, _out = _run(tmp_path, _prism_inp(with_cload=False, prescribed_top=True), "disp.inp")
    _uh, extra = results[0]
    sigma = extra["stress_cells"]
    for c in range(3):
        assert sigma[c, 2, 2] == pytest.approx(SIGMA_ZZ, abs=ABS_S)
        assert extra["principals"][c, 2] == pytest.approx(SIGMA_ZZ, abs=ABS_S)
        assert extra["von_mises"][c] == pytest.approx(abs(SIGMA_ZZ), abs=ABS_S)


def test_parse_only_nao_gera_manifesto_de_sucesso(tmp_path):
    payload, out = _run(tmp_path, _prism_inp(), "parse.inp", parse_only=True)
    assert payload["parse_only"] is True
    assert payload["report"]["status"] == "parse_ok_not_solved"
    assert payload["report"]["solver_run"] is False
    assert (out / "parse_report.json").is_file()
    status = json.loads((out / "run_status.json").read_text())
    assert status["status"] == "parse_ok_not_solved"
    assert status["solver_run"] is False
    assert not (out / "run_manifest.json").exists()


def test_solver_singular_nao_e_sucesso(tmp_path):
    # Sem BC, carga resultante nao nula: sistema singular / nao confinavel.
    text = _prism_inp(with_bc=False, with_cload=True)
    inp = tmp_path / "sing.inp"
    inp.write_text(text, encoding="utf-8")
    out = tmp_path / "out_sing"
    out.mkdir()
    stale = out / "lixo.pvd"
    stale.write_text("VELHO", encoding="utf-8")
    (out / "run_manifest.json").write_text(
        json.dumps({"status": "ok", "run_id": "old", "artifacts": [{"path": str(stale), "field": "U"}]})
        + "\n",
        encoding="utf-8",
    )
    with pytest.raises((RuntimeError, Exception)):
        run_from_inp(str(inp), output_dir=str(out), debug=0, profile="fc", units="mm-N-MPa-t-s")
    # A API direta nao grava run_status failed (isso e no __main__);
    # garante que nao reescreveu manifesto como sucesso desta run.
    if (out / "run_manifest.json").exists():
        man = json.loads((out / "run_manifest.json").read_text())
        assert man.get("run_id") in (None, "old") or man.get("status") != "ok" or man.get("n_elements") is None
        # manifesto de sucesso com n_elements seria desta run
        assert "n_elements" not in man or man.get("run_id") == "old"
    status_path = out / "run_status.json"
    if status_path.exists():
        st = json.loads(status_path.read_text())
        assert st["status"] in {"started", "failed"}


def test_cli_singular_exit_e_status_failed(tmp_path):
    import subprocess
    text = _prism_inp(with_bc=False)
    inp = tmp_path / "sing2.inp"
    inp.write_text(text, encoding="utf-8")
    out = tmp_path / "out_cli"
    out.mkdir()
    (out / "run_manifest.json").write_text('{"status": "ok", "run_id": "old"}\n', encoding="utf-8")
    src = os.path.join(ROOT, "src", "run_from_inp.py")
    proc = subprocess.run(
        [sys.executable, src, str(inp), "--profile", "fc", "--units", "mm-N-MPa-t-s", "--output-dir", str(out)],
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    st = json.loads((out / "run_status.json").read_text())
    assert st["status"] == "failed"
    assert st["run_id"] != "old"
    # manifesto antigo pode permanecer, mas status atual e failed
    man = json.loads((out / "run_manifest.json").read_text())
    assert man.get("run_id") == "old" or man.get("status") != "ok"



def test_axial_cg_gamg_vs_lu(tmp_path):
    """CG+GAMG no prisma axial deve ficar proximo do LU; relata o max |dU| e |dS|."""
    text = _prism_inp(split_cload=True)
    results_lu, out_lu = _run(tmp_path, text, "axial_lu.inp", solver="lu")
    results_cg, out_cg = _run(tmp_path, text, "axial_cg.inp", solver="cg-gamg", rtol=1e-10)
    _res_lu, extra_lu = results_lu[0]
    _res_cg, extra_cg = results_cg[0]
    import numpy as np
    u_lu = np.asarray(_res_lu.x.array, dtype=float)
    u_cg = np.asarray(_res_cg.x.array, dtype=float)
    du = float(np.max(np.abs(u_cg - u_lu)))
    s_lu = np.asarray(extra_lu["stress_cells"], dtype=float)
    s_cg = np.asarray(extra_cg["stress_cells"], dtype=float)
    ds = float(np.max(np.abs(s_cg - s_lu)))
    info = extra_cg.get("solver_info") or {}
    assert info.get("type") == "cg-gamg"
    assert int(info.get("converged_reason", 0)) > 0
    assert info.get("true_residual_rel") is not None
    assert info["true_residual_rel"] <= 1e-10
    assert info.get("true_residual_meets_rtol") is True
    # Nao inventa tolerancia de aceite para U/S vs LU; so relata a diferenca.
    print("AXIAL CG vs LU: max|dU|={0:.6e} max|dS|={1:.6e} its={2} ||r||/||b||={3:.6e}".format(
        du, ds, info.get("iterations"), info.get("true_residual_rel")))
    # consistencia grosseira: o prisma e minusculo; se CG funcionou, dU nao pode ser da ordem da solucao
    u_scale = float(np.max(np.abs(u_lu))) or 1.0
    assert du / u_scale < 1e-3
