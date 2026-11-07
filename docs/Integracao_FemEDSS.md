# Integração com FemEDSS

## Contrato atual (07/11/2026)

O FemEDSS-app isola cada análise em

```
uploads/<CNPJ>/<id do usuário>/<id do job>/{input,output}
```

e compacta só a pasta `output` daquele job. O FeniCS-EDSS precisa receber
o `.inp` e o diretório de saída por chamada. O contrato está em
`docs/Contrato_Integracao_FemEDSS.md`.

```bash
docker run --rm \
  -v "${FENICS_EDSS}:/app" \
  -v "${JOB}/input:/job/input:ro" \
  -v "${JOB}/output:/job/output" \
  edss/fenicsx:dev \
  python3 /app/src/run_from_inp.py /job/input/caso.inp --output-dir /job/output
```

Chamadas antigas, sem `--output-dir`, continuam gravando em `/app/output`.
O app não deve depender disso: dois jobs simultâneos iriam se misturar.

## Concorrência

`docker-compose.yaml` não define `container_name`. Vários `docker compose run`
ou `docker run` podem coexistir. Um compose de máquina com nome fixo de
container impede essa concorrência; o app deve passar um container novo por
job, ou usar Compose sem nome fixo.

O volume `./output:/app/output` permanece só para o padrão antigo. A saída
de um job do app é o volume montado em `/job/output`, não essa pasta
compartilhada.

## Compose chamado de dentro do container do app

Quando o `docker compose` roda via socket, a partir do container `app` do
FemEDSS, um volume relativo `.:/app` pode não apontar para o diretório certo
do host. Nesse caso o compose da máquina precisa montar o caminho absoluto
do FeniCS-EDSS no host. Esse arquivo local não é versionado aqui. Os
diretórios do job também precisam ser caminhos do host, visíveis no daemon
Docker, montados em `/job/input` e `/job/output`.
