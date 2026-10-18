# Templates Parametrizados — FeniCS-EDSS

## Visão geral

Os templates permitem executar simulações com parâmetros configuráveis via arquivo JSON, sem editar o código.

## Como usar

```bash
# Parâmetros padrão
./scripts/run_template.sh linear_elasticity

# Com arquivo de parâmetros
./scripts/run_template.sh linear_elasticity params/elasticity_default.json
```

## Templates disponíveis

### linear_elasticity
Elasticidade linear 3D (viga engastada sob peso próprio).

| Parâmetro | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| L | float | 1.0 | Comprimento da viga |
| W | float | 0.2 | Largura/altura da seção |
| nx, ny, nz | int | 20, 6, 6 | Divisões da malha |
| E | float | 1e5 | Módulo de Young |
| nu | float | 0.3 | Coeficiente de Poisson |
| rho | float | 1.0 | Densidade |
| g_scale | float | 0.4 | Escala da gravidade |
| output | string | /app/output/elasticity_param | Prefixo do arquivo de saída |

### heat
Equação do calor estacionária em domínio unitário.

| Parâmetro | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| nx, ny | int | 32, 32 | Divisões da malha |
| source | float | 10.0 | Intensidade da fonte de calor |
| output | string | /app/output/heat_param | Prefixo do arquivo de saída |

### hyperelasticity
Hiperelasticidade Neo-Hookean (viga sob carga).

| Parâmetro | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| L | float | 1.0 | Comprimento |
| W | float | 0.1 | Largura/altura |
| nx, ny, nz | int | 16, 4, 4 | Divisões da malha |
| E | float | 1e4 | Módulo de Young |
| nu | float | 0.3 | Coeficiente de Poisson |
| load_z | float | -0.5 | Carga na direção z |
| output | string | /app/output/hyper_param | Prefixo do arquivo de saída |

## Saída

Todos os templates geram:
- `.xdmf` + `.h5` (formato nativo DOLFINx)
- `.pvd` + `.vtu` (ParaView)

## Validação

Os templates validam:
- Tipos dos parâmetros
- Valores positivos onde aplicável (E, dimensões, divisões de malha)
- Intervalo de Poisson (0 < nu < 0.5)
