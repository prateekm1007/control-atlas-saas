FROM python:3.11-slim
RUN apt-get update && apt-get install -y curl && rm -rf /var/lib/apt/lists/*
WORKDIR /app
ENV PYTHONPATH=/app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY . /app
# Ensure all sub-packages are recognized
RUN find tos -type d -exec touch {}/__init__.py \;
EXPOSE 8000
# ENTRYPOINT BACK TO ROOT main.py
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
