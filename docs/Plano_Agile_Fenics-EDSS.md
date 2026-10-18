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

## Sprint 0 — Fundação (Concluído)
- Repositório GitHub criado
- Dockerfile baseado em `dolfinx/dolfinx:stable`
- docker-compose.yaml
- Exemplo Poisson funcional
- Saída em formato ParaView

## Sprint 1 — Solver mínimo (Concluído)
- Exemplos: Poisson, elasticidade linear, calor, hiperelasticidade
- Todos gerando `.pvd` estável no ParaView
- Script `scripts/run_fenics.sh`
- Commit e push realizados

## Sprint 2 — Templates parametrizados (Em andamento)
- Sistema de parâmetros via JSON
- Templates: `linear_elasticity`, `heat`, `hyperelasticity`
- Script `scripts/run_template.sh`
- Validação de parâmetros e logging
- Documentação em `docs/templates.md`

## Sprint 3 — Estabilização e testes (Planejado)
**Objetivo:** Garantir reprodutibilidade e qualidade.

- Testes automatizados (pytest) dos templates principais
- Validação de entrada mais robusta
- Tratamento de falhas (malha inválida, divergência, timeout)
- Ambiente de teste isolado no Nitro
- Versionamento de imagens Docker (tags `dev`, `test`)
- Pipeline simples de verificação local

**Critério de Done:** Todos os templates passam em testes automatizados no ambiente test do Nitro.

## Sprint 4 — Preparação para produção (Planejado)
**Objetivo:** FeniCS-EDSS rodando de forma confiável na Avell.

- Ajustar Dockerfile para produção (imagem menor, se possível)
- Definir estratégia de deploy na Avell (compose ou script)
- Testar execução com volumes e permissões corretas
- Definir limites de recursos adequados à Avell
- Documentar procedimento de atualização
- Tag de imagem `prod`

**Critério de Done:** Imagem de produção executa os templates principais na Avell sem intervenção manual.

## Sprint 5 — Integração com FemEDSS-app (Planejado)
**Objetivo:** FEniCSx disponível como terceira opção de solver no FemEDSS.

**Quando iniciar:** Após o Sprint 4 estar estável.

- Adicionar `fenicsx` como opção no modelo `Job` e na API
- Criar `executar_fenicsx()` em `app/services/solver.py`
- Ajustar interface web (seleção de solver + parâmetros)
- Garantir saída compatível (`.pvd` / `.vtu`)
- Testes de integração ponta a ponta
- Consumo de tokens funcionando

**Critério de Done:** Usuário do FemEDSS consegue submeter um job FEniCSx e baixar o resultado em ParaView.

## Princípios
- Docker desde o início
- Nitro = dev/test, Avell = prod
- Documentação sempre em `docs/`
- Commits frequentes na branch `dev`
