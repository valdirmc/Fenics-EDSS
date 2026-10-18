#!/usr/bin/env bash
# FeniCS-EDSS — executor de exemplos
# Uso: ./scripts/run_fenics.sh <exemplo>

set -e

EXEMPLO=${1:-}

if [ -z "$EXEMPLO" ]; then
    echo "Uso: $0 <exemplo>"
    echo ""
    echo "Exemplos disponíveis:"
    ls examples/*.py 2>/dev/null | xargs -n1 basename | sed 's/.py$//' || echo "  (nenhum)"
    exit 1
fi

SCRIPT="examples/${EXEMPLO}.py"

if [ ! -f "$SCRIPT" ]; then
    echo "Erro: arquivo não encontrado: $SCRIPT"
    exit 1
fi

echo "========================================"
echo " FeniCS-EDSS — $EXEMPLO"
echo "========================================"

mkdir -p output
rm -f output/${EXEMPLO}.*

echo "Executando..."
docker compose run --rm -T fenicsx python3 "/app/$SCRIPT"

echo ""
echo "Arquivos gerados:"
ls -lh output/${EXEMPLO}* 2>/dev/null || echo "(nenhum)"
