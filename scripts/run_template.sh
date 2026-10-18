#!/usr/bin/env bash
# Uso: ./scripts/run_template.sh <template> [arquivo_params.json]

set -e

TEMPLATE=${1:-}
PARAMS=${2:-}

if [ -z "$TEMPLATE" ]; then
    echo "Uso: $0 <template> [params.json]"
    echo "Templates disponíveis:"
    ls templates/*.py 2>/dev/null | xargs -n1 basename | sed 's/.py$//'
    exit 1
fi

SCRIPT="templates/${TEMPLATE}.py"
if [ ! -f "$SCRIPT" ]; then
    echo "Template não encontrado: $SCRIPT"
    exit 1
fi

echo "========================================"
echo " Template: $TEMPLATE"
echo "========================================"

if [ -n "$PARAMS" ] && [ -f "$PARAMS" ]; then
    echo "Usando parâmetros: $PARAMS"
    docker compose run --rm -T fenicsx python3 "/app/$SCRIPT" "/app/$PARAMS"
else
    echo "Usando parâmetros padrão"
    docker compose run --rm -T fenicsx python3 "/app/$SCRIPT"
fi