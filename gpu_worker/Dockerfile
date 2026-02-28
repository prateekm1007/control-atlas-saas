FROM nvidia/cuda:12.1.1-runtime-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# System dependencies
RUN apt-get update && apt-get install -y \
    python3.11 \
    python3-pip \
    python3.11-dev \
    wget \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Python dependencies
RUN pip3 install --no-cache-dir \
    openmm>=8.2 \
    pdbfixer>=1.12.0 \
    celery==5.3.6 \
    redis==5.0.1 \
    requests==2.31.0 \
    "numpy>=1.26.0,<3" 

WORKDIR /app

COPY worker/ /app/worker/
COPY protocols/ /app/protocols/

CMD ["celery", "-A", "worker.tasks", "worker", "--loglevel=info", "--concurrency=1"]
