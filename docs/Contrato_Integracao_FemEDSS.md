# Contrato de Integração FeniCS-EDSS ↔ FemEDSS-app

**Versão:** 1.0  
**Data:** 18/10/2026

## Visão geral

O FeniCS-EDSS é executado como container Docker. O FemEDSS-app envia um arquivo `.inp` e recebe os resultados em `.pvd`.

## Estrutura de pastas

```
Fenics-EDSS/
├── input/          # Arquivos .inp de entrada
├── output/         # Resultados .pvd gerados
├── src/            # Código do solver
└── docker-compose.yaml
```

## Chamada

```bash
# Dentro de Fenics-EDSS/
docker compose run --rm fenicsx python3 /app/src/run_from_inp.py /app/input/caso.inp [-debug=1]
```

## Entrada

- Arquivo `.inp` no formato CalculiX
- Colocado em `input/`
- Suporta: estática linear, térmica, hiperelasticidade, múltiplos steps e materiais

## Saída

- Arquivos `{nome}_{N}.pvd` em `output/` (um por step)
- Arquivos auxiliares `.pvtu`, `.vtu`

## Exemplo

```bash
cp meu_caso.inp input/
docker compose run --rm fenicsx python3 /app/src/run_from_inp.py /app/input/meu_caso.inp
ls output/meu_caso_*.pvd
```

## Integração com FemEDSS-app

O FemEDSS-app deve:
1. Salvar o `.inp` recebido em `Fenics-EDSS/input/`
2. Executar o comando Docker
3. Coletar os `.pvd` de `output/`
4. Disponibilizar para download ou visualização

## Volumes

O `docker-compose.yaml` já mapeia:
- `.:/app` (código e pastas input/output)
- `./output:/app/output`
