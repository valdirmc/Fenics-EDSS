# Plano Agile — FeniCS-EDSS (revisado)

**Projeto:** FeniCS-EDSS  
**Premissa:** O mesmo arquivo `.inp` que roda no CalculiX deve rodar no FeniCS-EDSS.  
**Data:** 18/10/2026

## Visão geral

| Sprint | Objetivo | Status |
|--------|----------|--------|
| 0–2 | Fundação, exemplos e templates básicos | Concluído |
| 3 | Parser `.inp` + Estática Linear | Em andamento |
| 4 | Condução Térmica + Não Linear | Planejado |
| 5 | Integração com FemEDSS-app | Planejado |

## Premissa central

A entrada principal é um arquivo `.inp` (formato CalculiX/Abaqus). O FeniCS-EDSS lê esse arquivo, traduz para DOLFINx e executa a análise.

## Cards `.inp` suportados na V1

### Malha
| Card | Descrição | V1 |
|------|-----------|----|
| `*NODE` | Coordenadas dos nós | Sim |
| `*ELEMENT` | Conectividade (C3D8, C3D4, C3D10) | Sim |

### Materiais e seções
| Card | Descrição | V1 |
|------|-----------|----|
| `*MATERIAL` | Definição de material | Sim |
| `*ELASTIC` | Isotrópico linear (E, nu) | Sim |
| `*SOLID SECTION` | Associação material → elementos | Sim |
| `*CONDUCTIVITY` | Condutividade térmica isotrópica | Sim (térmico) |
| `*HYPERELASTIC` | Neo-Hookean (C10, D1) | Sim (não linear) |

### Condições de contorno e cargas
| Card | Descrição | V1 |
|------|-----------|----|
| `*BOUNDARY` | Dirichlet (deslocamento ou temperatura) | Sim |
| `*CLOAD` | Força concentrada em nó | Sim |
| `*DLOAD` | Pressão/tração em face | Sim |
| `*DSLOAD` | Fonte de calor volumétrica | Sim (térmico) |

### Análise
| Card | Descrição | V1 |
|------|-----------|----|
| `*STATIC` | Estática linear ou incremental | Sim |
| `*HEAT TRANSFER` | Condução térmica estacionária | Sim |
| `*NLGEOM` | Ativa não linearidade geométrica | Sim |

### Controles (V1 simplificado)
| Card / Parâmetro | Descrição | V1 |
|------------------|-----------|----|
| Substeps | Número de incrementos | Sim (via parâmetro) |
| Save steps | Frequência de gravação | Sim (via parâmetro) |
| Solver type | Direto (LU) ou iterativo | Direto na V1 |

## Tipos de análise na V1

### 1. Estática Linear
- Materiais: `*ELASTIC` (múltiplos via `*SOLID SECTION`)
- Cargas: `*CLOAD`, `*DLOAD`
- BCs: `*BOUNDARY`
- Solver: direto

### 2. Condução Térmica
- Materiais: `*CONDUCTIVITY` (múltiplos)
- Fontes: `*DSLOAD` ou temperatura prescrita
- BCs: `*BOUNDARY` (temperatura)
- Solver: direto

### 3. Não Linear
- Materiais: `*ELASTIC` ou `*HYPERELASTIC`
- Geometria: `*NLGEOM`
- Cargas: múltiplos casos (incrementais)
- Controles: substeps, save steps
- Solver: Newton (direto no sistema linear)

## Elementos suportados na V1

| Tipo CalculiX | Célula DOLFINx | Dimensão |
|---------------|----------------|----------|
| C3D8 | Hexaedro | 3D |
| C3D4 | Tetraedro | 3D |
| C3D10 | Tetraedro quadrático | 3D |
| CPS4 / CPE4 | Quadrilátero | 2D (plane stress/strain) |

## Arquitetura

```
.inp
  │
  ▼
Parser (mesh, materials, BCs, loads, analysis type)
  │
  ▼
Tradutor → DOLFINx Mesh + FunctionSpace + Forms
  │
  ▼
Solver (LinearProblem ou NonlinearProblem)
  │
  ▼
.pvd / .vtu
```

## Sprint 3 — Parser + Estática Linear

- Parser mínimo de `.inp` (NODE, ELEMENT, MATERIAL, ELASTIC, SOLID SECTION, BOUNDARY, CLOAD, STATIC)
- Tradução para malha DOLFINx
- Estática linear com um material
- Saída `.pvd`

## Sprint 4 — Térmico + Não Linear

- Suporte a `*CONDUCTIVITY` e `*HEAT TRANSFER`
- Suporte a `*HYPERELASTIC` + `*NLGEOM`
- Múltiplos materiais
- Controles de substeps e save steps

## Sprint 5 — Integração FemEDSS

- `executar_fenicsx()` recebe caminho do `.inp`
- Retorna `.pvd`
- Compatível com o fluxo de jobs existente

## Princípios

- Mesmo `.inp` do CalculiX
- Docker
- Documentação em `docs/`
- Nitro = dev/test, Avell = prod
