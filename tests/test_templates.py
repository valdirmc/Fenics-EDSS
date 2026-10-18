"""
Testes básicos para os templates FeniCS-EDSS.
Rodar dentro do container: docker compose run --rm fenicsx python3 -m pytest /app/tests/ -v
"""

import json
import os
import sys
import pytest

# Adiciona o diretório raiz ao path
sys.path.insert(0, "/app")

def test_linear_elasticity_params_default():
    """Verifica que o template de elasticidade carrega parâmetros padrão."""
    from templates.linear_elasticity import load_params, validate
    p = load_params()
    assert p["L"] == 1.0
    assert p["E"] == 1.0e5
    validate(p)

def test_linear_elasticity_params_invalid_E():
    """Verifica que E negativo é rejeitado."""
    from templates.linear_elasticity import load_params, validate
    p = load_params()
    p["E"] = -100
    with pytest.raises(ValueError):
        validate(p)

def test_linear_elasticity_params_invalid_nu():
    """Verifica que Poisson fora do intervalo é rejeitado."""
    from templates.linear_elasticity import load_params, validate
    p = load_params()
    p["nu"] = 0.6
    with pytest.raises(ValueError):
        validate(p)

def test_heat_params_default():
    """Verifica parâmetros padrão do calor."""
    from templates.heat import load_params, validate
    p = load_params()
    assert p["nx"] == 32
    assert p["source"] == 10.0
    validate(p)

def test_heat_params_invalid_mesh():
    """Verifica que malha inválida é rejeitada."""
    from templates.heat import load_params, validate
    p = load_params()
    p["nx"] = 0
    with pytest.raises(ValueError):
        validate(p)

def test_hyperelasticity_params_default():
    """Verifica parâmetros padrão da hiperelasticidade."""
    from templates.hyperelasticity import load_params, validate
    p = load_params()
    assert p["L"] == 1.0
    assert p["load_z"] == -0.5
    validate(p)

def test_hyperelasticity_params_invalid_nu():
    """Verifica que Poisson inválido é rejeitado."""
    from templates.hyperelasticity import load_params, validate
    p = load_params()
    p["nu"] = -0.1
    with pytest.raises(ValueError):
        validate(p)

def test_params_json_exists():
    """Verifica que os arquivos de parâmetros de exemplo existem."""
    base = "/app/params"
    for name in ["elasticity_default.json", "heat_default.json", "hyper_default.json"]:
        path = os.path.join(base, name)
        # Aceita também o typo antigo
        if not os.path.exists(path):
            alt = path.replace("heat_default", "heat_deafult")
            if os.path.exists(alt):
                path = alt
        assert os.path.exists(path), f"Arquivo não encontrado: {path}"
        with open(path) as f:
            data = json.load(f)
            assert isinstance(data, dict)
