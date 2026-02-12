FROM python:3.11-bookworm
WORKDIR /app
RUN pip install --no-cache-dir streamlit requests
COPY . .
CMD ["streamlit", "run", "dashboard.py", "--server.port=8501", "--server.address", "0.0.0.0"]
