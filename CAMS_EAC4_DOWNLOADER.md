# CAMS EAC4 Downloader (Europa)

Script Python dockerizzato per scaricare CAMS global reanalysis (EAC4)
per area Europa, variabile per variabile e mese per mese.

## File

- Script: cams_eac4_download_europe.py
- Config: .env
- Compose service: cams-eac4-downloader

## Cosa fa

- Scarica un file per ogni (profilo, variabile, mese)
- Scrive subito su disco ogni chunk
- Supporta resume (skip file gia' validi)
- Retry con gestione throttling
- Report CSV finale

## Temporalita'

EAC4 e' 3-hourly:

- 00:00
- 03:00
- 06:00
- 09:00
- 12:00
- 15:00
- 18:00
- 21:00

## Profili supportati

- single: variabili single-level
- model: variabili multi-level su model levels
- pressure: variabili multi-level su pressure levels
- all: combina tutti i profili sopra

## Configurazione minima (.env)

Usa il blocco CAMS_* presente in .env.example.

Valori tipici richiesti:

- CAMS_DATASET=cams-global-reanalysis-eac4
- CAMS_START_YEAR=2003
- CAMS_END_YEAR=2025
- CAMS_AREA=72,-25,34,45
- CAMS_OUTPUT_DIR=/path/output
- CAMS_VARIABLES=all
- CAMS_LEVEL_TYPE=single
- CAMS_FILE_FORMAT=netcdf
- CAMSAPI_URL=https://ads.atmosphere.copernicus.eu/api
- CAMSAPI_KEY=<api-key-ads>

Nota formato:

- CAMS_FILE_FORMAT=grib -> file .grib
- CAMS_FILE_FORMAT=netcdf oppure netcdf_zip -> file .zip (netCDF zippato)

## Esecuzione Docker Compose

Dalla cartella satellites:

```bash
docker compose up cams-eac4-downloader
```

## Esecuzione docker run

```bash
docker build -t cams-downloader .
docker run --rm \
  --env-file .env \
  -e CAMS_OUTPUT_DIR=/data/cams_eac4 \
  -v "$(pwd)/data:/data" \
  --entrypoint python \
  cams-downloader cams_eac4_download_europe.py --env .env
```

## Note pratiche

- CAMS_VARIABLES=all usa il catalogo ADS per ottenere variabili disponibili nel profilo scelto.
- Se CAMS_LEVEL_TYPE=model o pressure e non imposti CAMS_MODEL_LEVELS/CAMS_PRESSURE_LEVELS,
  lo script richiede automaticamente tutti i livelli disponibili.
- Prima esecuzione consigliata: CAMS_DRY_RUN=true per verificare pianificazione e permessi.
