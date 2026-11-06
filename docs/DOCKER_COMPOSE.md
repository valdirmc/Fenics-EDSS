# Docker Compose — general vs por máquina

| Arquivo | No git? | Papel |
|---------|---------|-------|
| `docker-compose-general.yaml` | sim | Build, imagem `edss/fenicsx:dev`, volumes relativos `.:/app`, `sleep infinity`. |
| `docker-compose.yaml` | **não** (gitignored) | Paths absolutos (Nitro `/mnt/MAIN/...` ou `/run/media/...`; Avell `/home/valdir/...`) e overrides de nested Docker. |

## O que é geral vs específico

| Item | General | Por máquina |
|------|---------|-------------|
| build / Dockerfile / image | ✓ | |
| `command: sleep infinity` | ✓ | |
| volumes `.:/app`, `./output` | ✓ | |
| volume absoluto do repo | | ✓ (Nitro/Avell paths) |
| `FENICS_HOST_DIR` bind extra | | ✓ (Avell / nested Docker) |
| CPU/RAM limits | placeholder comentado | ✓ se precisar |

## Uso

```bash
# Overlay (recomendado)
docker compose -f docker-compose-general.yaml -f docker-compose.yaml up -d --build

# Ou copie o general e edite só o local:
cp docker-compose-general.yaml docker-compose.yaml
# edite volumes absolutos se o compose for chamado via docker.sock
docker compose up -d --build
```

## Risco no próximo pull (Avell)

O commit remove `docker-compose.yaml` do git. O próximo `git pull` na Avell **apaga** o arquivo trackeado do disco.
**Antes:** `cp docker-compose.yaml docker-compose.yaml.bak`  
**Depois:** restaure como `docker-compose.yaml` (agora ignorado) com o path `/home/valdir/...` se ainda precisar do bind absoluto.
