# FeniCS-EDSS — Documentação

## Visão geral

Projeto para integrar o FEniCSx como solver FEM no ecossistema EDSS.


## Documentos

| Documento | Conteúdo |
|-----------|----------|
| [ARCHITECTURE.md](ARCHITECTURE.md) | Overview técnico: fluxo `.inp`→DOLFINx→`.pvd`, EDSS vs terceiros, maturidade, bugs/riscos demo |
| [DOCKER_COMPOSE.md](DOCKER_COMPOSE.md) | Compose general vs por máquina |
| [Memoria_Projeto_Fenics-EDSS.md](Memoria_Projeto_Fenics-EDSS.md) | Memória cronológica / sprints |
| [Contrato_Integracao_FemEDSS.md](Contrato_Integracao_FemEDSS.md) | Contrato de integração com FemEDSS-app |
| [Integracao_FemEDSS.md](Integracao_FemEDSS.md) | Notas de volume/Docker aninhado |
| [Plano_Agile_Fenics-EDSS.md](Plano_Agile_Fenics-EDSS.md) | Plano Agile |
| [Cards_INP_V1.md](Cards_INP_V1.md) | Cards `.inp` V1 |
| [templates.md](templates.md) | Templates parametrizados |

## Ambientes

| Máquina | Função | Caminho |
|---------|--------|---------|
| Nitro | dev / test | `/run/media/valdir/MAIN/000_VALDIR/GIT_PJR/FemEDSS/Fenics-EDSS` |
| Avell | prod | `/home/valdir/000_VALDIR/00_CAE/FemEDSS/Fenics-EDSS` |

## Próximos passos (Sprint 0)

1. Criar repositório no GitHub
2. Clonar na pasta indicada no Nitro
3. Build da imagem Docker
4. Executar o exemplo Poisson
