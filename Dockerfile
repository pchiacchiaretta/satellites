FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /app

RUN apt-get update \
    && apt-get install -y --no-install-recommends libeccodes-tools \
    && rm -rf /var/lib/apt/lists/*

# Dipendenze minime per download ERA5 via cdsapi
COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r /app/requirements.txt

COPY . /app

ENTRYPOINT ["python", "era5_download_europe_hourly.py"]
CMD ["--env", ".env"]
