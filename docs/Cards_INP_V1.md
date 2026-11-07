# Cards `.inp` suportados na V1 — FeniCS-EDSS

**Data:** 18/10/2026  
**Status:** Implementado

## Malha
| Card | Suporte |
|------|---------|
| `*NODE` | Sim |
| `*ELEMENT` (C3D8, C3D4) | Sim |

## Materiais e seções
| Card | Suporte |
|------|---------|
| `*MATERIAL` | Sim |
| `*ELASTIC` | Sim |
| `*SOLID SECTION` | Sim |
| `*ELSET` | Sim |
| `*NSET` | Sim |

## Condições de contorno e cargas
| Card | Suporte |
|------|---------|
| `*BOUNDARY` (nó ou nset, range de DOFs) | Sim |
| `*CLOAD` | Sim |
| `*DLOAD` | Parser ok, aplicação pendente |

## Análise
| Card | Suporte |
|------|---------|
| `*STATIC` | Sim |
| `*NLGEOM` | Sim |
| `*HEAT TRANSFER` | Sim |
| `*HYPERELASTIC` (Neo-Hooke: mu, lambda) | Sim |
| `*CONDUCTIVITY` | Sim |

## Status da implementação (18/10/2026)

| Recurso | Parser | Solver |
|---------|--------|--------|
| `*NSET` + `*BOUNDARY` | OK | OK |
| `*CLOAD` | OK | OK |
| `*DLOAD` | OK | OK (aproximação) |
| Múltiplos materiais | OK | OK |
| `*HEAT TRANSFER` | OK | OK |
| `*NLGEOM` / Hiperelasticidade | OK | OK (Neo-Hooke) |

## Exemplos
- `examples/simple_beam.inp` — estática, BCs por nó
- `examples/beam_nset.inp` — estática, BCs por nset
- `examples/beam_multimat.inp` — 2 materiais (parser)
- `examples/beam_dload.inp` — pressão
- `examples/beam_thermal.inp` — condução térmica (`*HEAT TRANSFER`)
- `examples/beam_hyper.inp` — hiperelasticidade (`*HYPERELASTIC`)

## Uso
```bash
# Saída no diretório do job
python3 /app/src/run_from_inp.py /job/input/beam_thermal.inp --output-dir /job/output

# Compatível com a chamada antiga (grava em /app/output)
python3 /app/src/run_from_inp.py /app/examples/beam_thermal.inp -debug=1
```
