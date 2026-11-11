"""Perfil FC: whitelist, rejeicoes e parse de cartoes (sem DOLFINx, sem roda)."""

import os
import sys

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


