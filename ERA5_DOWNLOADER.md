# ERA5 Hourly Single Levels — Downloader per Europa

Script Python dockerizzato per scaricare dati di rianalisi oraria ERA5
dal servizio **Copernicus Climate Data Store (CDS)** per l'intera area europea,
variabile per variabile e mese per mese.

---

## Indice

1. [Cosa è ERA5](#1-cosa-è-era5)
2. [Cosa fa lo script](#2-cosa-fa-lo-script)
3. [Struttura dei file](#3-struttura-dei-file)
4. [Prerequisiti](#4-prerequisiti)
5. [Configurazione (.env)](#5-configurazione-env)
6. [Variabili disponibili](#6-variabili-disponibili)
7. [Prima di lanciare](#7-prima-di-lanciare)
8. [Lancio con Docker Compose](#8-lancio-con-docker-compose-consigliato)
9. [Lancio con docker run](#9-lancio-con-docker-run)
10. [Struttura output](#10-struttura-output)
11. [Log — come leggere e monitorare](#11-log--come-leggere-e-monitorare)
12. [Report CSV riassuntivo](#12-report-csv-riassuntivo)
13. [Re-lancio con periodo più ampio](#13-re-lancio-con-periodo-più-ampio)
14. [Controlli automatici](#14-controlli-automatici)
15. [Troubleshooting](#15-troubleshooting)

---

## 1. Cosa è ERA5

**ERA5** è il dataset di rianalisi climatica globale prodotto da **ECMWF**
(European Centre for Medium-Range Weather Forecasts) nell'ambito del programma
**Copernicus Climate Change Service (C3S)**.

Caratteristiche principali:

| Proprietà | Valore |
|-----------|--------|
| Copertura temporale | 1940 – presente |
| Risoluzione spaziale | ~31 km (0.25°) |
| Risoluzione temporale | Oraria |
| Distribuzione | Copernicus CDS (accesso libero con registrazione) |
| Formati | GRIB2, NetCDF |

I dati vengono distribuiti tramite il **CDS API** (`cdsapi`), che accoda
le richieste e le evade in background; lo script gestisce automaticamente
attesa, retry e resume.

---

## 2. Cosa fa lo script

```
era5_download_europe_hourly.py
```

Per ogni **variabile** × ogni **mese** del periodo configurato:

1. **Valida** la richiesta prima di inviarla al CDS (parametri, area, formato, path).
2. **Controlla** se il file è già presente su disco (verifica anche l'header GRIB/NetCDF).
   - Se il file è valido → lo **salta** senza fare nessuna richiesta CDS.
   - Se il file è corrotto o vuoto → lo **elimina e riscarica**.
3. **Invia** la richiesta al CDS e attende il completamento.
4. **Verifica** il file scaricato (esistenza, dimensione, magic bytes).
5. **Scrive immediatamente** il file nella directory di output.
6. **Aggiorna** il log e il report CSV con lo stato del chunk.

Il download è **incrementale**: ogni file viene salvato subito, indipendentemente
dagli altri, senza aspettare il completamento dell'intero job.

---

## 3. Struttura dei file

```
satellites/
├── era5_download_europe_hourly.py   # Script principale
├── .env                             # Configurazione (NON committare la API key)
├── Dockerfile                       # Immagine Docker
├── docker-compose.yml               # Compose per lancio semplice
├── .dockerignore                    # Esclude file inutili dall'immagine
├── requirements.txt                 # Dipendenze Python (incluso cdsapi)
└── ERA5_DOWNLOADER.md               # Questa documentazione
```

Output generato a runtime:

```
<ERA5_OUTPUT_DIR>/
├── 2m_temperature/
│   ├── era5_2m_temperature_201401.grib
│   ├── era5_2m_temperature_201402.grib
│   └── ...
├── total_precipitation/
│   └── ...
├── era5_downloader.log              # Log completo di ogni run
└── era5_report.csv                  # Report CSV riassuntivo
```

---

## 4. Prerequisiti

### Account CDS

1. Registrati su [https://cds.climate.copernicus.eu](https://cds.climate.copernicus.eu).
2. Accetta i **Terms of Use** del dataset
   `reanalysis-era5-single-levels` (una volta sola, dalla pagina del dataset).
3. Vai sul tuo profilo → sezione **API key** → copia la chiave.

### Software host

| Requisito | Versione minima |
|-----------|-----------------|
| Docker    | 20.x            |
| Docker Compose | v2 (`docker compose`) |

Nessuna dipendenza Python è richiesta sull'host.

---

## 5. Configurazione (.env)

Il file `.env` nella cartella `satellites/` contiene tutti i parametri.
Editalo prima di lanciare lo script.

```dotenv
# ── DATASET ───────────────────────────────────────────────
ERA5_DATASET=reanalysis-era5-single-levels
ERA5_PRODUCT_TYPE=reanalysis

# ── PERIODO (inclusivo) ───────────────────────────────────
ERA5_START_YEAR=2014
ERA5_START_MONTH=1
ERA5_END_YEAR=2026
ERA5_END_MONTH=12

# ── AREA GEOGRAFICA [nord, ovest, sud, est] ───────────────
# Europa completa:
ERA5_AREA=72,-25,34,45

# ── OUTPUT ────────────────────────────────────────────────
ERA5_OUTPUT_DIR=/home/manager/storage/atmo/.../era5

# ── FORMATO FILE ─────────────────────────────────────────
ERA5_FILE_FORMAT=grib          # grib oppure netcdf

# ── VARIABILI (nomi CDS, separati da virgola) ────────────
ERA5_VARIABLES=2m_temperature,total_precipitation,...

# ── RETRY E CONTROLLO FLUSSO ─────────────────────────────
ERA5_MAX_RETRIES=3
ERA5_RETRY_WAIT_SECONDS=25
ERA5_REQUEST_PAUSE_SECONDS=2
ERA5_SKIP_EXISTING=true
ERA5_DRY_RUN=false

# ── LOGGING ───────────────────────────────────────────────
ERA5_LOG_FILE=                 # vuoto = <ERA5_OUTPUT_DIR>/era5_downloader.log
ERA5_LOG_LEVEL=INFO            # DEBUG, INFO, WARNING, ERROR
ERA5_VERIFY_DOWNLOAD_MAGIC=true

# ── REPORT CSV ────────────────────────────────────────────
ERA5_REPORT_FILE=              # vuoto = <ERA5_OUTPUT_DIR>/era5_report.csv

# ── CREDENZIALI CDS ───────────────────────────────────────
CDSAPI_URL=https://cds.climate.copernicus.eu/api
CDSAPI_KEY=<la-tua-api-key>
```

### Tabella parametri

| Parametro | Descrizione | Default |
|-----------|-------------|---------|
| `ERA5_DATASET` | Nome dataset CDS | `reanalysis-era5-single-levels` |
| `ERA5_PRODUCT_TYPE` | Tipo prodotto | `reanalysis` |
| `ERA5_START_YEAR` / `ERA5_END_YEAR` | Anno inizio/fine (inclusivi) | `2014` / `2026` |
| `ERA5_START_MONTH` / `ERA5_END_MONTH` | Mese inizio/fine (1–12) | `1` / `12` |
| `ERA5_AREA` | Bounding box `nord,ovest,sud,est` | `72,-25,34,45` |
| `ERA5_OUTPUT_DIR` | Cartella destinazione file | **obbligatorio** |
| `ERA5_FILE_FORMAT` | `grib` oppure `netcdf` | `grib` |
| `ERA5_VARIABLES` | Lista variabili CDS (virgola) | 16 variabili default |
| `ERA5_MAX_RETRIES` | Tentativi per chunk in caso di errore | `3` |
| `ERA5_RETRY_WAIT_SECONDS` | Attesa tra un retry e il successivo (s) | `25` |
| `ERA5_REQUEST_PAUSE_SECONDS` | Pausa tra richieste successive al CDS (s) | `2` |
| `ERA5_SKIP_EXISTING` | Salta file già presenti e validi | `true` |
| `ERA5_DRY_RUN` | Simula senza scaricare nulla | `false` |
| `ERA5_LOG_FILE` | Percorso log (vuoto = nella output dir) | — |
| `ERA5_LOG_LEVEL` | Livello log | `INFO` |
| `ERA5_VERIFY_DOWNLOAD_MAGIC` | Verifica header GRIB/NetCDF post-download | `true` |
| `ERA5_REPORT_FILE` | Percorso CSV riassuntivo (vuoto = nella output dir) | — |
| `CDSAPI_URL` | URL API CDS | `https://cds.climate.copernicus.eu/api` |
| `CDSAPI_KEY` | API key del tuo account CDS | **obbligatorio** |

---

## 6. Variabili disponibili

Le 16 variabili configurate di default:

| Nome CDS | Descrizione |
|----------|-------------|
| `2m_temperature` | Temperatura a 2 m |
| `10m_u_component_of_wind` | Componente U del vento a 10 m |
| `10m_v_component_of_wind` | Componente V del vento a 10 m |
| `convective_precipitation` | Precipitazione convettiva |
| `convective_rain_rate` | Tasso di pioggia convettiva |
| `instantaneous_10m_wind_gust` | Raffica di vento istantanea a 10 m |
| `k_index` | Indice K (instabilità atmosferica) |
| `precipitation_type` | Tipo di precipitazione |
| `total_cloud_cover` | Copertura nuvolosa totale |
| `total_column_ozone` | Ozono colonna totale |
| `total_column_rain_water` | Acqua di pioggia in colonna totale |
| `total_column_snow_water` | Acqua nevosa in colonna totale |
| `total_column_supercooled_liquid_water` | Acqua liquida sopraffusa in colonna totale |
| `total_column_water` | Acqua totale in colonna |
| `total_column_water_vapour` | Vapore acqueo in colonna totale |
| `total_precipitation` | Precipitazione totale |

Per aggiungere o rimuovere variabili modifica `ERA5_VARIABLES` nel `.env`.
L'elenco completo dei nomi CDS è disponibile su:
[https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels)

---

## 7. Prima di lanciare

### 1. Imposta la API key nel .env

```
CDSAPI_KEY=<la-tua-api-key>
```

### 2. Imposta il percorso di output

```
ERA5_OUTPUT_DIR=/percorso/dove/salvare/i/file
```

Il percorso **deve essere raggiungibile dal container** (vedi sezione Compose).

### 3. Verifica lo spazio disponibile

Un mese di dati ERA5 orari per 1 variabile su Europa pesa tipicamente
**2–8 GB** in formato GRIB. Con 16 variabili × 156 mesi (2014–2026) prevedi
**nell'ordine di qualche TB** totali.

Controlla lo spazio libero:

```bash
df -h /percorso/dove/salvare/i/file
```

### 4. (Opzionale) Dry run per testare la configurazione

```dotenv
ERA5_DRY_RUN=true
```

Esegue tutto il ciclo di validazione e log senza inviare nessuna richiesta
al CDS e senza scrivere file. Utile per verificare che i parametri siano
corretti. Rimetti `false` prima del lancio reale.

---

## 8. Lancio con Docker Compose (consigliato)

### Prima volta / rebuild immagine

```bash
cd /home/manager/services/satellites
docker compose build
```

### Avvio

```bash
docker compose up
```

Per mandarlo in background:

```bash
docker compose up -d
```

### Seguire i log in tempo reale

```bash
# log del container (stdout)
docker compose logs -f era5-downloader

# oppure dal file di log (più comodo per sessioni lunghe)
tail -f /home/manager/storage/atmo/sftpgo/inbound/raw_data/satellites/copernicus/era5/era5_downloader.log
```

### Verificare che l'area sia davvero Europa

Controlla i limiti geografici direttamente nel file GRIB (non solo nel `.env`):

```bash
docker compose exec era5-downloader sh -lc \
'grib_ls -w count=1 -p latitudeOfFirstGridPointInDegrees,longitudeOfFirstGridPointInDegrees,latitudeOfLastGridPointInDegrees,longitudeOfLastGridPointInDegrees /data/era5/2m_temperature/era5_2m_temperature_201401.grib'
```

Con `ERA5_AREA=72,-25,34,45` l'output corretto deve mostrare:

- `latitudeOfFirstGridPointInDegrees = 72`
- `longitudeOfFirstGridPointInDegrees = -25`
- `latitudeOfLastGridPointInDegrees = 34`
- `longitudeOfLastGridPointInDegrees = 45`

Se questi 4 valori coincidono, stai scaricando esattamente il dominio Europa impostato.

### Stop

```bash
docker compose down
```

Il download in corso viene interrotto pulitamente; al prossimo `up` riprende
dal primo file non ancora scaricato (grazie a `ERA5_SKIP_EXISTING=true`).

---

### Note sul volume di output

Nel `docker-compose.yml` il percorso di output interno al container è
`/data/era5`, montato sulla cartella host reale:
`/home/manager/storage/atmo/sftpgo/inbound/raw_data/satellites/copernicus/era5`.

```yaml
volumes:
  - /home/manager/storage/atmo/sftpgo/inbound/raw_data/satellites/copernicus/era5:/data/era5
environment:
  ERA5_OUTPUT_DIR: /data/era5
```

Se vuoi scrivere su un percorso host diverso (es. un disco esterno),
modifica il lato sinistro del bind mount:

```yaml
volumes:
  - /mnt/disco_esterno/era5:/data/era5
```

---

## 9. Lancio con docker run

```bash
cd /home/manager/services/satellites

docker build -t era5-downloader .

docker run --rm \
  --name era5-downloader \
  --env-file .env \
  -e ERA5_OUTPUT_DIR=/data/era5 \
  -v "/home/manager/storage/atmo/sftpgo/inbound/raw_data/satellites/copernicus/era5:/data/era5" \
  era5-downloader --env .env
```

---

## 10. Struttura output

I file vengono scritti **immediatamente** al completamento di ogni chunk,
organizzati per variabile:

```
<ERA5_OUTPUT_DIR>/
├── 2m_temperature/
│   ├── era5_2m_temperature_201401.grib
│   ├── era5_2m_temperature_201402.grib
│   ├── ...
│   └── era5_2m_temperature_202612.grib
├── 10m_u_component_of_wind/
│   └── ...
├── total_precipitation/
│   └── ...
├── era5_downloader.log
└── era5_report.csv
```

**Convenzione nome file:**

```
era5_<variabile>_<YYYY><MM>.<ext>
```

Esempi:
- `era5_2m_temperature_202403.grib`
- `era5_total_precipitation_201506.nc`

---

## 11. Log — come leggere e monitorare

Il log viene scritto sia su **stdout del container** sia su **file**
(`era5_downloader.log` nella cartella di output).

### Formato

```
2026-03-29T17:30:01 | INFO | [REQUEST_READY] dataset=reanalysis-era5-single-levels variable=2m_temperature year=2014 month=01 ...
2026-03-29T17:30:01 | INFO | [DOWNLOAD_START] variable=2m_temperature year=2014 month=01 attempt=1/3
2026-03-29T17:35:42 | INFO | [DOWNLOAD_DONE] variable=2m_temperature year=2014 month=01 target=... ok, size=47123456
2026-03-29T17:35:44 | INFO | [SKIP_EXISTING_OK] variable=2m_temperature year=2014 month=02 ...
```

### Tabella eventi

| Evento | Livello | Significato |
|--------|---------|-------------|
| `[REQUEST_READY]` | INFO | Richiesta validata, pronta per il CDS |
| `[DOWNLOAD_START]` | INFO | Download avviato (con numero tentativo) |
| `[DOWNLOAD_DONE]` | INFO | File scaricato, verificato e salvato |
| `[SKIP_EXISTING_OK]` | INFO | File già presente e valido, saltato |
| `[SKIP_EXISTING_BAD]` | WARNING | File presente ma corrotto → verrà riscaricato |
| `[DOWNLOAD_ERROR]` | ERROR | Errore durante il download (con dettaglio) |
| `[VERIFY_FAILED]` | ERROR | File scaricato ma header non valido |
| `[RETRY_WAIT]` | INFO | Attesa prima di riprovare |
| `[TARGET_INVALID]` | ERROR | Path o estensione del file non coerente |
| `[REQUEST_INVALID]` | ERROR | Parametri della richiesta non coerenti |
| `[DRY_RUN]` | INFO | Dry-run attivo, nessun download effettuato |
| `[REPORT]` | INFO | Report CSV scritto a fine sessione |

### Comandi utili per monitorare

```bash
# Segui il log in tempo reale
tail -f /percorso/output/era5_downloader.log

# Conta quanti file sono stati scaricati con successo
grep -c "DOWNLOAD_DONE" /percorso/output/era5_downloader.log

# Conta quanti file saltati (già presenti)
grep -c "SKIP_EXISTING_OK" /percorso/output/era5_downloader.log

# Vedi solo gli errori
grep "\[DOWNLOAD_ERROR\]\|\[VERIFY_FAILED\]\|\[TARGET_INVALID\]\|\[REQUEST_INVALID\]" \
  /percorso/output/era5_downloader.log

# Vedi i segnali di throttling/rate-limit
grep "\[THROTTLE_DETECTED\]\|\[THROTTLE_WAIT\]\|\[ADAPTIVE_PAUSE_UP\]\|\[ADAPTIVE_PAUSE_DOWN\]" \
  /percorso/output/era5_downloader.log

# Ultima variabile in lavorazione
grep "DOWNLOAD_START\|DOWNLOAD_DONE" /percorso/output/era5_downloader.log | tail -5
```

---

## 12. Report CSV riassuntivo

Al termine di ogni run lo script scrive (o aggiorna in **append**) un file CSV:

```
<ERA5_OUTPUT_DIR>/era5_report.csv
```

### Colonne

| Colonna | Descrizione |
|---------|-------------|
| `timestamp` | Data e ora del chunk (ISO 8601) |
| `variable` | Nome variabile ERA5 |
| `year` | Anno |
| `month` | Mese (2 cifre) |
| `status` | `downloaded` / `skipped` / `error` / `dry_run` |
| `throttle_detected` | `true` se nel chunk e' stato rilevato rate-limit/throttling |
| `file_path` | Percorso assoluto del file |
| `file_size_bytes` | Dimensione file in byte |
| `note` | Dettaglio (es. `ok, size=47123456` oppure messaggio di errore) |

### Esempio

```csv
timestamp,variable,year,month,status,throttle_detected,file_path,file_size_bytes,note
2026-03-29T17:30:01,2m_temperature,2014,01,downloaded,false,/data/.../era5_2m_temperature_201401.grib,47123456,ok size=47123456
2026-03-29T17:35:44,2m_temperature,2014,02,skipped,false,/data/.../era5_2m_temperature_201402.grib,46980000,ok size=46980000
2026-03-29T17:36:01,convective_rain_rate,2014,01,error,true,/data/.../era5_convective_rain_rate_201401.grib,0,max retries esauriti: ...
```

### Analisi rapida del report

```bash
# Riepilogo conteggi per status
awk -F',' 'NR>1 {print $5}' era5_report.csv | sort | uniq -c

# Lista file in errore
awk -F',' 'NR>1 && $5=="error" {print $2, $3, $4, $8}' era5_report.csv

# Totale GB scaricati
awk -F',' 'NR>1 && $5=="downloaded" {sum+=$7} END {printf "%.2f GB\n", sum/1024/1024/1024}' era5_report.csv
```

---

## 13. Re-lancio con periodo più ampio

Se in futuro vuoi estendere il periodo (es. da 2014–2026 a 2010–2027),
basta modificare il `.env` e rilanciare: i file già scaricati vengono
**riconosciuti e saltati** automaticamente.

```dotenv
ERA5_START_YEAR=2010   # era 2014
ERA5_END_YEAR=2027     # era 2026
```

```bash
docker compose up
```

Lo script:
1. Per ogni chunk controlla se il file esiste già nella cartella.
2. Se esiste **e l'header è valido** → evento `SKIP_EXISTING_OK`, nessuna chiamata al CDS.
3. Se esiste ma è corrotto → evento `SKIP_EXISTING_BAD`, lo elimina e riscarica.
4. Se non esiste → scarica normalmente.

> ⚠️ Il comportamento di skip dipende da `ERA5_SKIP_EXISTING=true` nel `.env`.
> Se lo metti a `false` tutti i file vengono riscaricati sovrascrivendo quelli esistenti.

---

## 14. Controlli automatici

Lo script esegue i seguenti controlli prima, durante e dopo ogni download:

### Pre-download

| Controllo | Cosa verifica |
|-----------|---------------|
| Path target | Il file è dentro `ERA5_OUTPUT_DIR` e ha l'estensione corretta |
| Payload request | `dataset`, `variable`, `year`, `month`, `day`, `time` (24 ore), `area`, `data_format` coerenti |

### Post-download

| Controllo | Cosa verifica |
|-----------|---------------|
| Esistenza file | Il file esiste su disco dopo il download |
| Dimensione | Il file non è vuoto |
| Header magic bytes | GRIB → inizia con `GRIB`; NetCDF → inizia con `CDF` o `\x89HDF` |

Se la verifica post-download fallisce, il file viene eliminato e
il download viene ritentato (fino a `ERA5_MAX_RETRIES` volte).

> Disabilitare la verifica: `ERA5_VERIFY_DOWNLOAD_MAGIC=false`

---

## 15. Troubleshooting

### Errore: `missing API key` o `401 Unauthorized`

Verifica che `CDSAPI_KEY` nel `.env` sia valorizzata senza spazi:

```dotenv
CDSAPI_KEY=99cf1230-e27f-4620-b330-cfc6e0dc08b7   # corretto
CDSAPI_KEY= 99cf1230-e27f-4620-b330-cfc6e0dc08b7  # sbagliato (spazio iniziale)
```

### Errore: `Terms of Use not accepted`

Accedi su [https://cds.climate.copernicus.eu](https://cds.climate.copernicus.eu),
vai alla pagina del dataset `reanalysis-era5-single-levels` e accetta i
Terms of Use (operazione una tantum).

### Il download va in timeout o errore di rete

La coda CDS può essere molto affollata. Lo script riprova automaticamente
(`ERA5_MAX_RETRIES`). Se i retry si esauriscono, il chunk viene marcato
`error` nel CSV e il log mostra `[DOWNLOAD_ERROR]`. Al lancio successivo
quei chunk verranno ritentati (il file non è stato salvato, quindi non
vengono saltati).

### Un file nel CSV è marcato `error` ma voglio forzare il retry

Basta rilanciare lo script: se il file non è presente su disco (o è vuoto),
viene automaticamente ritentato.

### Spazio disco esaurito durante il download

Il file parziale (vuoto o incompleto) viene eliminato dallo script.
Libera spazio e rilancia: i file già scaricati correttamente vengono saltati.

### Verificare quanti file mancano senza scaricarli

Imposta `ERA5_DRY_RUN=true` e rilancia: il log mostrerà un `[DRY_RUN]`
per ogni chunk che verrebbe scaricato, e `[SKIP_EXISTING_OK]` per quelli
già presenti. Poi rimetti `ERA5_DRY_RUN=false`.
