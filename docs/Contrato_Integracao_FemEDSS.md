# Contrato de Integração FeniCS-EDSS ↔ FemEDSS-app

**Versão:** 1.1
**Data:** 07/11/2026

## Visão geral

O FeniCS-EDSS roda como container Docker, um por job. O FemEDSS-app entrega
o `.inp` daquele job e lê somente a pasta de saída daquele job. O app compacta
essa pasta para download.

Layout no FemEDSS-app:

```
uploads/<CNPJ>/<id do usuário>/<id do job>/input/<arquivo>.inp
uploads/<CNPJ>/<id do usuário>/<id do job>/output/
```

## Linha de comando

```bash
python3 /app/src/run_from_inp.py <arquivo.inp> [--output-dir DIR] [-debug=1]
```

| Argumento | Obrigatório | Significado |
|-----------|-------------|-------------|
| `arquivo.inp` | sim | Caminho do `.inp` deste job, dentro do container |
| `--output-dir DIR` | não | Diretório de saída deste job. Atalhos: `-output DIR`, `-output=DIR`, `--output-dir=DIR` |
| `-debug=1` | não | Diagnóstico no stdout. Padrão: `0` |

Se `--output-dir` não for informado, a saída vai para `/app/output`.
Esse padrão existe para os scripts antigos. O FemEDSS-app **deve** passar
o diretório do job. Uma opção desconhecida encerra com erro, em vez de
gravar na pasta global.

## O que o app deve executar

Subir um container novo por job. Não usar `container_name` fixo e não
fazer `docker exec` num container compartilhado.

```bash
JOB="/caminho/uploads/<CNPJ>/<id do usuário>/<id do job>"

docker run --rm \
  -v "${FENICS_EDSS}:/app" \
  -v "${JOB}/input:/job/input:ro" \
  -v "${JOB}/output:/job/output" \
  edss/fenicsx:dev \
  python3 /app/src/run_from_inp.py \
    "/job/input/<arquivo>.inp" \
    --output-dir /job/output
```

`${FENICS_EDSS}` é o diretório deste repositório no host (a imagem não copia
`src/`). `${JOB}/output` precisa existir ou ser criável pelo Docker antes
do processo gravar; o solver também cria esse diretório se o caminho
montado permitir.

Com Compose, o equivalente é `docker compose run --rm -T` **sem**
`container_name`, acrescentando os dois `-v` do job. O
`docker-compose.yaml` versionado não define nome fixo de container.
Um compose local de máquina, com caminho absoluto do host, não faz parte
deste contrato e não deve ser exigido pelo app.

## Saída

Somente dentro de `--output-dir`:

- `{nome do inp}_{número do step}.pvd`
- `{nome}_{número}000000.pvtu` e `{nome}_{número}_p0_000000.vtu` (nomes gerados pelo DOLFINx)
- na condução térmica, também `.xdmf` e `.h5` ao lado do `.pvd`

O `.pvd` aponta para o `.pvtu` com caminho relativo ao próprio arquivo, e o
`.pvtu` aponta para o `.vtu` do mesmo jeito. O zip de `output/` abre no
ParaView sem a árvore do servidor.

O solver não apaga arquivos, não faz glob em `/app/output` nem em outras
pastas de job, e não reescreve `.pvd` que ele não acabou de gerar.

## Entrada `.inp`

Formato CalculiX/Abaqus já aceito pela V1: estática linear, `*HEAT TRANSFER`
(condução), `*HYPERELASTIC` (Neo-Hooke, com ou sem `*NLGEOM`), vários steps
e materiais. O tipo do step escolhe o solver:

- `*HEAT TRANSFER` → térmico
- `*HYPERELASTIC` ou `*NLGEOM` → hiperelástico
- caso contrário → estática linear

## Exemplo local (um job)

```bash
mkdir -p /tmp/job-demo/input /tmp/job-demo/output
cp examples/beam_nset.inp /tmp/job-demo/input/
docker run --rm \
  -v "$PWD:/app" \
  -v /tmp/job-demo/input:/job/input:ro \
  -v /tmp/job-demo/output:/job/output \
  edss/fenicsx:dev \
  python3 /app/src/run_from_inp.py /job/input/beam_nset.inp --output-dir /job/output
ls /tmp/job-demo/output
```
