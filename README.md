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
├── docker-compose-general.yaml  # versionado
├── docker-compose.yaml         # local (gitignored)
├── README.md
└── .gitignore
```

## Ambientes

- **Nitro** (dev/test): `/run/media/valdir/MAIN/000_VALDIR/GIT_PJR/FemEDSS/Fenics-EDSS`
- **Avell** (prod): `/home/valdir/000_VALDIR/00_CAE/FemEDSS/Fenics-EDSS`


## Docker Compose

- Shared: [`docker-compose-general.yaml`](docker-compose-general.yaml)
- Por máquina: `docker-compose.yaml` (gitignored)
- Detalhes: [docs/DOCKER_COMPOSE.md](docs/DOCKER_COMPOSE.md)

```bash
docker compose -f docker-compose-general.yaml -f docker-compose.yaml build
docker compose -f docker-compose-general.yaml -f docker-compose.yaml run --rm fenicsx python3 /app/examples/poisson.py
```
## Uso rápido

```bash
# Build
docker compose build

# Executar exemplo Poisson
docker compose run --rm fenicsx python3 /app/examples/poisson.py
```

## Documentação

- Overview técnico (fluxo, maturidade, riscos): [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md)
- Índice: [docs/README.md](docs/README.md)
