"""Perfil FC: whitelist, rejeicoes e parse de cartoes (sem DOLFINx, sem roda)."""

import os
import sys
from pathlib import Path

import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))

import numpy as np
from fc_profile import FCProfileError, normalize_units  # noqa: E402
from fc_fields import principals_algebraic  # noqa: E402
from inp_parser import parse_inp  # noqa: E402
from job_io import parse_run_args  # noqa: E402


def _w(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return str(path)


MIN_C3D4 = """*NODE, NSET=Nall
1, 0, 0, 0
2, 1, 0, 0
3, 0, 1, 0
4, 0, 0, 1
*ELEMENT, TYPE=C3D4, ELSET=Evolumes
1, 1, 2, 3, 4
*ELSET, ELSET=Eall
Evolumes
*NSET, NSET=Fixos
1, 2, 3
*MATERIAL, NAME=Mat
*ELASTIC
70000, 0.33
*SOLID SECTION, ELSET=Eall, MATERIAL=Mat
*STEP
*STATIC
1.0, 1.0
*BOUNDARY
Fixos,1
Fixos,2
Fixos,3
*CLOAD
** comentario no meio do cartao
4, 3, -1.0
*NODE FILE
U
*EL FILE
S, E
*NODE PRINT, NSET=Fixos, TOTALS=ONLY
RF
*END STEP
"""


def test_fc_aceita_c3d4_estatico_e_bc_duas_colunas(tmp_path):
    model = parse_inp(_w(tmp_path, "ok.inp", MIN_C3D4), profile="fc")
    assert len(model.elements) == 1
    assert next(iter(model.elements.values())).type == "C3D4"
    assert len(model.nsets["Fixos"]) == 3
    assert model.elsets["Eall"] == [1]
    assert {bc.dof for bc in model.boundary} == {1, 2, 3}
    assert len(model.boundary) == 9
    assert model.steps[0].loads[0].value == -1.0
    fields = []
    for req in model.fc_record.requests:
        fields.extend(req.fields)
    assert "U" in fields and "S" in fields and "E" in fields and "RF" in fields
    assert model.steps[0].start_time == 0.0
    assert model.steps[0].static_increments == [1.0, 1.0]


def test_fc_rejeita_type_ausente(tmp_path):
    text = MIN_C3D4.replace("TYPE=C3D4, ", "")
    with pytest.raises(FCProfileError, match="sem TYPE"):
        parse_inp(_w(tmp_path, "no_type.inp", text), profile="fc")


def test_fc_rejeita_c3d8(tmp_path):
    text = MIN_C3D4.replace("TYPE=C3D4", "TYPE=C3D8").replace("1, 1, 2, 3, 4", "1, 1, 2, 3, 4, 1, 2, 3, 4")
    with pytest.raises(FCProfileError, match="C3D4"):
        parse_inp(_w(tmp_path, "hex.inp", text), profile="fc")


def test_fc_rejeita_multistep(tmp_path):
    text = MIN_C3D4.replace("*END STEP", "*END STEP\n*STEP\n*STATIC\n*END STEP")
    with pytest.raises(FCProfileError, match="um \\*STEP"):
        parse_inp(_w(tmp_path, "multi.inp", text), profile="fc")


@pytest.mark.parametrize(
    "card,frag",
    [
        ("*DLOAD\nEall, P1, 1.0\n", "DLOAD"),
        ("*NLGEOM\n", "NLGEOM"),
        ("*HEAT TRANSFER\n", "HEAT TRANSFER"),
        ("*HYPERELASTIC, NEO HOOKE\n1.0, 0.1\n", "HYPERELASTIC"),
        ("*PLASTIC\n", "PLASTIC"),
        ("*FREQUENCY\n", "FREQUENCY"),
        ("*TEMPERATURE\n1, 20\n", "TEMPERATURE"),
        ("*INITIAL CONDITIONS, TYPE=TEMPERATURE\nNall, 20\n", "INITIAL CONDITIONS"),
        ("*FOO\n", "whitelist"),
    ],
)
def test_fc_rejeita_bloqueados_e_desconhecidos(tmp_path, card, frag):
    text = MIN_C3D4.replace("*END STEP", card + "*END STEP")
    # *HYPERELASTIC precisa de material corrente; coloca depois de *ELASTIC
    if "HYPERELASTIC" in card:
        text = MIN_C3D4.replace("*SOLID SECTION", card + "*SOLID SECTION")
    with pytest.raises(FCProfileError, match=frag):
        parse_inp(_w(tmp_path, "bad.inp", text), profile="fc")


def test_fc_ignora_metadados_explicitos(tmp_path):
    text = MIN_C3D4.replace(
        "*SOLID SECTION",
        "*PHYSICAL CONSTANTS, ABSOLUTE ZERO=0\n*EXPANSION, ZERO=0\n2.36E-05\n*SOLID SECTION",
    ).replace("*STEP", "*STEP, INC=2000")
    model = parse_inp(_w(tmp_path, "meta.inp", text), profile="fc")
    keys = {item.keyword for item in model.fc_record.ignored}
    assert "PHYSICAL CONSTANTS" in keys
    assert "EXPANSION" in keys
    assert any(item.keyword == "STEP" for item in model.fc_record.ignored)


def test_perfil_default_ainda_aceita_heat_e_type_default(tmp_path):
    model = parse_inp(_w(tmp_path, "def.inp", "*ELEMENT\n1, 1, 2, 3, 4, 5, 6, 7, 8\n*STEP\n*HEAT TRANSFER\n*END STEP\n"))
    assert next(iter(model.elements.values())).type == "C3D8"
    assert model.analysis_type == "heat"


def test_cli_profile_e_units():
    args = parse_run_args(["caso.inp", "--profile", "fc", "--units", "mm-N-s-t", "--output-dir", "/tmp/out"])
    assert args.profile == "fc"
    assert args.units == "mm-N-s-t"
    assert normalize_units(args.units) == "mm-N-MPa-t-s"
    with pytest.raises(ValueError, match="perfil desconhecido"):
        parse_run_args(["caso.inp", "--profile", "xyz"])


def test_principais_algebricos_com_sinal():
    # Tensor com tracao e compressao.
    sigma = np.zeros((1, 3, 3))
    sigma[0, 0, 0] = -10.0
    sigma[0, 1, 1] = 4.0
    sigma[0, 2, 2] = 1.0
    pr, vm = principals_algebraic(sigma)
    assert pr[0, 0] == pytest.approx(4.0)
    assert pr[0, 1] == pytest.approx(1.0)
    assert pr[0, 2] == pytest.approx(-10.0)
    assert vm[0] == pytest.approx((0.5 * ((4-1)**2 + (1-(-10))**2 + ((-10)-4)**2)) ** 0.5)


def test_fc_rejeita_output_frequency_nao_honrada(tmp_path):
    text = MIN_C3D4.replace("*NODE FILE", "*OUTPUT, FREQUENCY=2\n*NODE FILE")
    with pytest.raises(FCProfileError, match="OUTPUT nao honrado"):
        parse_inp(_w(tmp_path, "out2.inp", text), profile="fc")


def test_fc_aceita_output_frequency_1(tmp_path):
    text = MIN_C3D4.replace("*NODE FILE", "*OUTPUT, FREQUENCY=1\n*NODE FILE")
    model = parse_inp(_w(tmp_path, "out1.inp", text), profile="fc")
    keys = [k for _ln, k in ((item[0], item[1]) for item in model.fc_record.keywords)]
    # keywords store header strings
    headers = [h for _i, h in model.fc_record.keywords]
    assert any("OUTPUT" in h.upper() for h in headers)


def test_fc_rejeita_expansion_com_temperatura(tmp_path):
    text = MIN_C3D4.replace(
        "*SOLID SECTION",
        "*EXPANSION, ZERO=0\n2.36E-05\n*TEMPERATURE\n1, 20\n*SOLID SECTION",
    )
    with pytest.raises(FCProfileError, match="TEMPERATURE"):
        parse_inp(_w(tmp_path, "exp_temp.inp", text), profile="fc")


def test_parse_nao_afirma_solver_valido(tmp_path):
    from fc_profile import PARSE_NOT_SOLVED_BANNER, parse_report_payload
    model = parse_inp(_w(tmp_path, "ok.inp", MIN_C3D4), profile="fc")
    payload = parse_report_payload(model, "ok.inp")
    assert payload["status"] == "parse_ok_not_solved"
    assert payload["solver_run"] is False
    assert "somente de parse" in payload["warning"]
    assert "parse_ok_not_solved" in PARSE_NOT_SOLVED_BANNER
    assert "Nao implica solver valido" in PARSE_NOT_SOLVED_BANNER


def test_cli_parse_only():
    args = parse_run_args(["caso.inp", "--profile", "fc", "--units", "mm-N-MPa-t-s", "--parse-only", "--output-dir", "/tmp/out"])
    assert args.parse_only is True


def test_run_status_e_lixo_velho_nao_contam(tmp_path):
    from job_io import write_run_manifest, write_run_status, verify_run_artifacts, file_sha256
    out = tmp_path / "out"
    out.mkdir()
    stale = out / "old.pvd"
    stale.write_text("VELHO", encoding="utf-8")
    old_manifest = write_run_manifest(str(out), {
        "status": "ok",
        "run_id": "old-run",
        "artifacts": [{"path": str(stale), "field": "U"}],
    })
    # nova run falha: status failed, mesmo manifesto antigo no disco
    status = write_run_status(str(out), {"status": "failed", "run_id": "new-run", "error": "teste"})
    data = (out / "run_status.json").read_text(encoding="utf-8")
    assert '"status": "failed"' in data
    assert "new-run" in data
    # manifesto antigo continua no disco, mas o status atual nao e ok
    assert Path(old_manifest).is_file()
    # hash do lixo velho nao deve ser usado como sucesso da nova run
    stale.write_text("ALTERADO", encoding="utf-8")
    with pytest.raises(RuntimeError, match="lixo velho"):
        verify_run_artifacts(old_manifest, ["old.pvd"])


def test_principais_compressao_uniaxial_e_cisalhamento():
    sigma = np.zeros((2, 3, 3))
    # celula 0: compressao uniaxial
    sigma[0, 2, 2] = -100.0
    # celula 1: cisalhamento puro + residual assimetrico
    sigma[1, 0, 1] = 30.0
    sigma[1, 1, 0] = 29.0
    pr, vm = principals_algebraic(sigma)
    assert pr[0, 0] == pytest.approx(0.0)
    assert pr[0, 1] == pytest.approx(0.0)
    assert pr[0, 2] == pytest.approx(-100.0)
    assert vm[0] == pytest.approx(100.0)
    assert pr[1, 0] >= pr[1, 1] >= pr[1, 2]
    # cisalhamento simetrizado (~29.5): principais +29.5, 0, -29.5
    assert pr[1, 0] == pytest.approx(29.5)
    assert pr[1, 1] == pytest.approx(0.0)
    assert pr[1, 2] == pytest.approx(-29.5)
    assert vm[1] == pytest.approx((0.5 * (2 * (29.5 ** 2) + (29.5 - (-29.5)) ** 2)) ** 0.5)


def test_fc_rejeita_boundary_3_colunas_valor(tmp_path):
    # Formato "no, dof, valor" e ambiguo: CalculiX interpreta a 3a coluna como dof_fim.
    text = MIN_C3D4.replace("Fixos,1\nFixos,2\nFixos,3", "4, 3, -0.001")
    with pytest.raises(FCProfileError, match="3 colunas"):
        parse_inp(_w(tmp_path, "bc3.inp", text), profile="fc")


def test_fc_aceita_boundary_4_colunas_valor(tmp_path):
    text = MIN_C3D4.replace("Fixos,1\nFixos,2\nFixos,3", "4, 3, 3, -0.001")
    model = parse_inp(_w(tmp_path, "bc4.inp", text), profile="fc")
    vals = {(bc.node, bc.dof, bc.value) for bc in model.boundary}
    assert (4, 3, -0.001) in vals


def test_t0_hidrostatico_von_mises():
    """C6/T0 proposto: tensao hidostatica deve ter sigma_vm ~ 0 (sem cancelamento)."""
    for pval in (100.0, -100.0, 1.0):
        sigma = np.eye(3)[None, :, :] * pval
        pr, vm = principals_algebraic(sigma)
        assert pr[0, 0] == pytest.approx(pval)
        assert pr[0, 2] == pytest.approx(pval)
        assert vm[0] <= 1e-12 * abs(pval)
