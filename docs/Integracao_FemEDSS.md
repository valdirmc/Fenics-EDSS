mkdir -p docs
cat > docs/Integracao_FemEDSS.md << 'EOF'
# Integração com FemEDSS

## Status (19/10/2026)
- Volume absoluto necessário quando o `docker compose` é chamado de dentro do container `app` do FemEDSS (Docker aninhado via socket).
- Volume relativo `.:/app` não resolve corretamente nesse cenário.

## Configuração atual (`docker-compose.yaml`)
```yaml
volumes:
  - /home/valdir/000_VALDIR/00_CAE/Fenics-EDSS:/app
  - /home/valdir/000_VALDIR/00_CAE/Fenics-EDSS/output:/app/output

 ## Chamada a partir do FemEDSS
  O FenicsService executa docker compose run --rm fenicsx python3 /app/src/run_from_inp.py ... com cwd em /fenics.
   Resultados gravados em output/.

text**2. Commit:**

```bash
git add docker-compose.yaml docs/Integracao_FemEDSS.md
git commit -m "fix: volume absoluto para integração com FemEDSS

- Substitui volume relativo por caminho absoluto do host
- Necessário quando o compose é chamado via docker.sock a partir do container app
- Documentação em docs/Integracao_FemEDSS.md"