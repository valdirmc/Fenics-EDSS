# FeniCS-EDSS

Solver FEM baseado em FEniCSx (DOLFINx) para o ecossistema EDSS / FemEDSS.

## Objetivo

Oferecer o FEniCSx como terceira opção de solver (ao lado de CalculiX e OpenCourant), com execução isolada em Docker.

## Estrutura

```
Fenics-EDSS/
├── docker/
│   └── Dockerfile
├── examples/
│   └── poisson.py
├── scripts/
│   └── run_fenics.sh
├── docs/
├── tests/
├── docker-compose.yaml
├── README.md
└── .gitignore
```

## Ambientes

- **Nitro** (dev/test): `/run/media/valdir/MAIN/000_VALDIR/GIT_PJR/FemEDSS/Fenics-EDSS`
- **Avell** (prod): `/home/valdir/000_VALDIR/00_CAE/FemEDSS/Fenics-EDSS`

## Uso rápido

```bash
# Build
docker compose build

# Executar exemplo Poisson
docker compose run --rm fenicsx python3 /app/examples/poisson.py
```
