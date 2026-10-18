# Cards `.inp` suportados na V1 — FeniCS-EDSS

**Data:** 18/10/2026

## Malha

| Card | Sintaxe mínima | Suporte V1 |
|------|----------------|------------|
| `*NODE` | `*NODE` + linhas `id, x, y, z` | Sim |
| `*ELEMENT` | `*ELEMENT, TYPE=C3D8` + conectividade | Sim (C3D8, C3D4, C3D10) |

## Materiais

| Card | Sintaxe mínima | Suporte V1 |
|------|----------------|------------|
| `*MATERIAL` | `*MATERIAL, NAME=nome` | Sim |
| `*ELASTIC` | `*ELASTIC` + `E, nu` | Sim |
| `*CONDUCTIVITY` | `*CONDUCTIVITY` + `k` | Sim (térmico) |
| `*HYPERELASTIC` | `*HYPERELASTIC, NEO HOOKE` + `C10, D1` | Sim (não linear) |
| `*SOLID SECTION` | `*SOLID SECTION, ELSET=..., MATERIAL=...` | Sim |

## Condições de contorno e cargas

| Card | Sintaxe mínima | Suporte V1 |
|------|----------------|------------|
| `*BOUNDARY` | `*BOUNDARY` + `nó, dof, valor` | Sim |
| `*CLOAD` | `*CLOAD` + `nó, dof, força` | Sim |
| `*DLOAD` | `*DLOAD` + `elset, P, pressão` | Sim |
| `*DSLOAD` | `*DSLOAD` + fonte volumétrica | Sim (térmico) |

## Análise

| Card | Sintaxe mínima | Suporte V1 |
|------|----------------|------------|
| `*STATIC` | `*STATIC` | Sim |
| `*HEAT TRANSFER` | `*HEAT TRANSFER, STEADY STATE` | Sim |
| `*NLGEOM` | `*NLGEOM` | Sim |

## Não suportados na V1

- `*DYNAMIC`
- `*FREQUENCY`
- `*PLASTIC`
- `*CONTACT`
- Elementos de casca, viga, treliça
- Análise acoplada termo-mecânica

## Exemplo mínimo (estática linear)

```
*NODE
1, 0.0, 0.0, 0.0
2, 1.0, 0.0, 0.0
...
*ELEMENT, TYPE=C3D8
1, 1, 2, 3, 4, 5, 6, 7, 8
...
*MATERIAL, NAME=Aco
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=Eall, MATERIAL=Aco
*BOUNDARY
1, 1, 1, 0.0
1, 2, 2, 0.0
1, 3, 3, 0.0
*CLOAD
10, 3, -100.0
*STEP
*STATIC
*END STEP
```
