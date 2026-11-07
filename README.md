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

## Solver a partir de .inp

Cada job informa o arquivo de entrada e o diretório de saída daquele job.
Se `--output-dir` for omitido, a saída continua em `/app/output`.

```bash
python3 /app/src/run_from_inp.py /job/input/caso.inp --output-dir /job/output
```

Formas aceitas da opção: `--output-dir DIR`, `--output-dir=DIR`, `-output DIR`, `-output=DIR`.
O diagnóstico segue `-debug=1`.

O solver grava `{nome}_{passo}.pvd`, `.pvtu` e `.vtu` somente nesse diretório
(na térmica, também `.xdmf`/`.h5`). As referências do `.pvd` para os `.vtu`
ficam relativas, para o zip da pasta abrir no ParaView. Não há limpeza nem
busca de arquivos em outras pastas.

### Jobs simultâneos

O `docker-compose.yaml` deste repositório não fixa `container_name`.
O FemEDSS-app deve subir um container por job, com o layout
`uploads/<CNPJ>/<id do usuário>/<id do job>/{input,output}`:

```bash
docker run --rm \
  -v "${FENICS_EDSS}:/app" \
  -v "${JOB}/input:/job/input:ro" \
  -v "${JOB}/output:/job/output" \
  edss/fenicsx:dev \
  python3 /app/src/run_from_inp.py /job/input/caso.inp --output-dir /job/output
```

Não reutilize um container de nome fixo (`docker exec fenicsx-edss`).
Dois jobs ao mesmo tempo precisam de diretórios de saída diferentes.
O contrato completo está em `docs/Contrato_Integracao_FemEDSS.md`.
