# Contrato de Integração — FeniCS-EDSS ↔ FemEDSS-app

**Data:** 18/10/2026  
**Status:** Proposta

## Objetivo

Definir como o FemEDSS-app vai chamar o FeniCS-EDSS como terceiro solver (ao lado de CalculiX e OpenCourant).

## Visão geral

```
FemEDSS-app (FastAPI)
    │
    ├── Job queue
    │
    └── executar_fenicsx()
            │
            └── docker compose / docker run
                    └── FeniCS-EDSS container
                            ├── Template (heat, elasticity, hyperelasticity)
                            ├── Parâmetros JSON
                            └── Saída: .pvd / .vtu
```

## Entrada

### Opção A — Template + parâmetros (recomendada)

O FemEDSS envia:

- Nome do template (`linear_elasticity`, `heat`, `hyperelasticity`)
- Arquivo JSON com parâmetros (ou dicionário)
- Diretório de trabalho (onde a saída será gravada)

Exemplo de chamada:

```python
executar_fenicsx(
    template="linear_elasticity",
    params={"E": 2.0e5, "nu": 0.3, "nx": 16},
    workdir="/caminho/do/job"
)
```

### Opção B — Script livre

O usuário envia um script Python completo. O FeniCS-EDSS apenas executa. Menos controlado, útil para casos avançados.

## Saída

| Item | Formato | Observação |
|------|---------|------------|
| Resultado principal | `.pvd` + `.vtu` | Compatível com ParaView |
| Resultado nativo | `.xdmf` + `.h5` | Formato DOLFINx |
| Status | string | `concluido` ou `erro` |
| Mensagem | string | Log resumido ou mensagem de erro |

## Contrato da função no FemEDSS-app

```python
def executar_fenicsx(
    template: str,
    params: dict,
    workdir: str,
    timeout: int = 600
) -> tuple[str, str]:
    """
    Executa um template FeniCS-EDSS.
    
    Returns:
        (status, caminho_resultado)
        status: "concluido" ou "erro"
        caminho_resultado: caminho do .pvd ou mensagem de erro
    """
```

## Implementação sugerida em `app/services/solver.py`

```python
def executar_fenicsx(template: str, params: dict, workdir: str, timeout: int = 600) -> tuple[str, str]:
    import json
    import os
    import subprocess
    
    # Grava parâmetros em arquivo temporário
    params_file = os.path.join(workdir, "params.json")
    with open(params_file, "w") as f:
        json.dump(params, f)
    
    # Define saída
    params["output"] = os.path.join(workdir, template)
    
    cmd = [
        "docker", "compose", "run", "--rm", "-T",
        "-v", f"{workdir}:/app/output",
        "-v", f"{params_file}:/app/params.json",
        "fenicsx",
        "python3", f"/app/templates/{template}.py", "/app/params.json"
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, cwd=workdir)
        pvd = os.path.join(workdir, f"{template}.pvd")
        if result.returncode == 0 and os.path.exists(pvd):
            return "concluido", pvd
        return "erro", result.stderr or result.stdout or "Falha na execução"
    except subprocess.TimeoutExpired:
        return "erro", "Timeout na execução"
    except Exception as e:
        return "erro", str(e)
```

## Ajustes necessários no FemEDSS-app

1. **Modelo `Job`**: adicionar `fenicsx` como valor possível em `solver`
2. **API**: aceitar `template` e `params` quando `solver == "fenicsx"`
3. **Interface web**: seletor de template + formulário de parâmetros (ou upload de JSON)
4. **Serviço**: implementar `executar_fenicsx()` conforme acima

## Volume e permissões

- O diretório do job no host deve ser montado em `/app/output` no container
- O container precisa de permissão de escrita nesse volume
- Na Avell, o caminho base será `/home/valdir/000_VALDIR/00_CAE/FemEDSS/`

## Próximos passos

1. Validar este contrato
2. Implementar `executar_fenicsx()` no FemEDSS-app (Sprint 5)
3. Testar ponta a ponta com um job real
