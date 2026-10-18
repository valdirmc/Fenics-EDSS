# Plano Agile — FeniCS-EDSS

**Projeto:** FeniCS-EDSS  
**Objetivo:** Solver FEM baseado em FEniCSx como 3ª opção no FemEDSS  
**Data:** 18/10/2026

## Visão geral dos Sprints

| Sprint | Objetivo | Ambiente | Status |
|--------|----------|----------|--------|
| 0 | Fundação do repositório e ambiente Docker | Nitro | Concluído |
| 1 | Solver mínimo (4 exemplos + script) | Nitro | Concluído |
| 2 | Templates parametrizados | Nitro | Em andamento |
| 3 | Estabilização e testes | Nitro test | Planejado |
| 4 | Preparação para produção | Avell | Planejado |
| 5 | Integração com FemEDSS-app | Nitro + Avell | Planejado |

## Sprint 0 — Fundação
- Repositório GitHub
- Dockerfile + docker-compose
- Exemplo Poisson funcional

## Sprint 1 — Solver mínimo
- Exemplos: Poisson, elasticidade linear, calor, hiperelasticidade
- Saída `.pvd` estável
- Script `run_fenics.sh`

## Sprint 2 — Templates parametrizados
- Sistema de parâmetros JSON
- Templates: linear_elasticity, heat, hyperelasticity
- Script `run_template.sh`
- Validação e logging (em melhoria)

## Sprint 5 — Integração (quando retomar o FemEDSS-app)
- Adicionar `fenicsx` como opção de solver no modelo Job
- Criar `executar_fenicsx()` em `app/services/solver.py`
- Manter saída compatível (`.pvd` / `.vtu`)

## Princípios
- Docker desde o início
- Nitro = dev/test, Avell = prod
- Documentação sempre em `docs/`
- Commits frequentes na branch `dev`
