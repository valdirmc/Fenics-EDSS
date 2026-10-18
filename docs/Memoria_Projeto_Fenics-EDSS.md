# Memória do Projeto FeniCS-EDSS

**Data de atualização:** 18/10/2026  
**Status:** Sprint 2 em andamento

## Visão geral

O **FeniCS-EDSS** é o terceiro solver FEM do ecossistema EDSS/FemEDSS, baseado em FEniCSx (DOLFINx). Roda isolado em Docker e será integrado ao FemEDSS-app como opção de solver (ao lado de CalculiX e OpenCourant).

## Ambientes

| Máquina | Função | Caminho |
|---------|--------|---------|
| Nitro | dev / test | `/run/media/valdir/MAIN/000_VALDIR/GIT_PJR/FemEDSS/Fenics-EDSS` |
| Avell | prod | `/home/valdir/000_VALDIR/00_CAE/FemEDSS/Fenics-EDSS` |

Repositório: https://github.com/valdirmc/Fenics-EDSS (branch `dev`)

## Histórico

### Sprint 0 — Fundação
- Estrutura inicial do repositório
- Dockerfile baseado em `dolfinx/dolfinx:stable`
- docker-compose.yaml
- Exemplo Poisson funcional
- Saída em formato ParaView

### Sprint 1 — Solver mínimo
- Exemplos: Poisson, elasticidade linear, calor, hiperelasticidade
- Todos gerando `.pvd` (estável no ParaView)
- Script `scripts/run_fenics.sh` para execução padronizada
- Commit e push realizados

### Sprint 2 — Templates parametrizados (em andamento)
- Templates criados: `linear_elasticity`, `heat`, `hyperelasticity`
- Sistema de parâmetros via JSON
- Script `scripts/run_template.sh`
- Parâmetros de exemplo em `params/`
- Validação de parâmetros e logging adicionados

### Sprint 3 — Parser `.inp` + Análises (concluído)
- Parser de `.inp` com suporte a:
  - `*NODE`, `*ELEMENT` (C3D8, C3D4)
  - `*MATERIAL`, `*ELASTIC`, `*CONDUCTIVITY`
  - `*NSET`, `*ELSET`, `*SOLID SECTION`
  - `*BOUNDARY` (nó ou nset, range de DOFs)
  - `*CLOAD`, `*DLOAD`
  - `*STATIC`, `*HEAT TRANSFER`, `*NLGEOM`
- Reordenação de hexaedros CalculiX → DOLFINx/Basix
- Estática linear funcional
- Condução térmica funcional
- Debug controlado por `-debug=1`
- Múltiplos materiais por elset (DG-0) funcional
- V1 completa: estática, térmica, hiperelasticidade, múltiplos steps e materiais

## Estrutura atual

```
Fenics-EDSS/
├── docker/
│   └── Dockerfile
├── examples/
├── templates/
├── params/
├── scripts/
├── docs/
├── output/
├── docker-compose.yaml
└── README.md
```

## Licença do FEniCSx

Componentes principais sob **LGPL-3.0-or-later** (DOLFINx, FFCx, UFL) e **MIT** (Basix). Compatível com uso comercial.

## Próximos passos

- Documentar templates
- Preparar integração com FemEDSS-app
- Deploy em produção na Avell
