# Base image with Python 3.10
FROM python:3.10-slim

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy project files
COPY . /app/control-atlas

# Install Python dependencies (PINNED for RDKit compatibility)
RUN pip install --no-cache-dir \
    "numpy<2" \
    rdkit-pypi \
    prody

# Set Python path explicitly (no undefined var usage)
ENV PYTHONPATH=/app/control-atlas

# Set entrypoint to Universal Gatekeeper
ENTRYPOINT ["python3", "/app/control-atlas/entries/020_candidate_validation/universal_gatekeeper.py"]

# Default command
CMD ["--help"]
