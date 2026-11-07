"""Cards .inp que precisam rotear: *HEAT TRANSFER e *HYPERELASTIC."""

import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))

from inp_parser import (  # noqa: E402
    get_material_for_element,
    parse_inp,
    select_analysis,
)

EXAMPLES = os.path.join(ROOT, "examples")


def _escreve(tmp_path, nome, texto):
    path = tmp_path / nome
    path.write_text(texto, encoding="utf-8")
    return str(path)


def test_heat_transfer_roteia_termico():
    model = parse_inp(os.path.join(EXAMPLES, "beam_thermal.inp"))
    assert model.analysis_type == "heat"
    assert model.steps
    assert model.steps[0].analysis_type == "heat"
    assert select_analysis(model, model.steps[0]) == "heat"
    materiais = {m.name.upper(): m for m in model.materials.values()}
    assert materiais["ACO"].k == 50.0
    assert materiais["ALUMINIO"].k == 200.0


def test_hyperelastic_roteia_e_le_parametros():
    model = parse_inp(os.path.join(EXAMPLES, "beam_hyper.inp"))
    assert model.hyperelastic is True
    assert model.nlgeom is True
    borracha = next(m for m in model.materials.values() if m.hyperelastic)
    assert borracha.mu == 1.0
    assert borracha.lmbda == 0.1
    assert model.steps
    assert select_analysis(model, model.steps[-1]) == "hyper"


def test_hiperelastico_sem_nlgeom_ainda_roteia(tmp_path):
    path = _escreve(
        tmp_path,
        "neo.inp",
        """*MATERIAL, NAME=Borracha
*HYPERELASTIC, NEO HOOKE
2.5, 0.2
*STEP
*STATIC
*END STEP
""",
    )
    model = parse_inp(path)
    assert model.hyperelastic is True
    assert model.nlgeom is False
    mat = next(iter(model.materials.values()))
    assert mat.mu == 2.5
    assert mat.lmbda == 0.2
    assert select_analysis(model, model.steps[0]) == "hyper"


def test_estatica_linear_continua_linear():
    model = parse_inp(os.path.join(EXAMPLES, "simple_beam.inp"))
    assert model.analysis_type == "static"
    assert model.hyperelastic is False
    for step in model.steps:
        assert select_analysis(model, step) == "linear"


def test_get_material_for_element_retorna_fallback_e_elset(tmp_path):
    path = _escreve(
        tmp_path,
        "mats.inp",
        """*ELEMENT, TYPE=C3D8
1, 1, 2, 3, 4, 5, 6, 7, 8
2, 1, 2, 3, 4, 5, 6, 7, 8
*ELSET, ELSET=Esquerda
1
*MATERIAL, NAME=Aco
*ELASTIC
210000, 0.3
*MATERIAL, NAME=Aluminio
*ELASTIC
70000, 0.33
*SOLID SECTION, ELSET=Esquerda, MATERIAL=Aluminio
""",
    )
    model = parse_inp(path)
    especifico = get_material_for_element(model, 1)
    assert especifico is not None
    assert especifico.name.upper() == "ALUMINIO"
    assert especifico.E == 70000
    fallback = get_material_for_element(model, 2)
    assert fallback is not None
    assert fallback.name.upper() == "ACO"
    vazio = parse_inp(_escreve(tmp_path, "vazio.inp", "*STEP\n*END STEP\n"))
    assert get_material_for_element(vazio, 1) is None
