FROM dolfinx/dolfinx:stable

LABEL maintainer="EDSS / Valdir"
LABEL description="FEniCSx solver environment for EDSS"

# Dependências extras úteis
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    vim \
    && rm -rf /var/lib/apt/lists/*

# Diretório de trabalho
WORKDIR /app

# Copia exemplos e scripts (serão montados em runtime na maioria dos casos)
COPY examples/ /app/examples/

# Variáveis úteis
ENV PYTHONUNBUFFERED=1
ENV PYVISTA_OFF_SCREEN=true

# Ponto de entrada padrão
CMD ["python3"]
