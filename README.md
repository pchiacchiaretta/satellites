# **CAMS – Mappe annuali, medie regionali, GIF, montage e time-series per capoluoghi**

**Abruzzo + Puglia**

Questo progetto genera:

* **mappe annuali** di inquinanti atmosferici (da CAMS)
* **medie regionali pesate**
* **GIF** multi-anno
* **montage multi-anno** (anni affiancati con **un’unica colorbar**)
* **grafici di andamento temporale (time-series) per capoluogo di provincia**

---

## 1) Cos’è CAMS e cosa stai usando qui

**CAMS (Copernicus Atmosphere Monitoring Service)** è il servizio Copernicus dedicato al monitoraggio atmosferico: produce e distribuisce **analisi**, **previsioni** e **reanalisi** (ricostruzioni storiche consistenti) per vari componenti atmosferici e per la qualità dell’aria.

Nel tuo caso utilizzi i file NetCDF della **reanalisi/produzione regionale europea** di qualità dell’aria (dataset mensili *validated* o *interim* e/o file già annuali).
Questi prodotti combinano modellistica e osservazioni e sono pensati per analisi climatologiche, trend temporali e confronti spaziali.

---

## 2) Cosa fa lo script (in breve)

Lo script:

1. **Scansiona i file NetCDF** in `DATA_DIR`

   * Se trova file annuali:

     ```
     annual_mean_<pollutant>_<year>.nc
     ```

     li usa direttamente (più veloce).
   * Altrimenti calcola la media annua dai mensili:

     ```
     cams-europe-air-quality-reanalyses_<pollutant>_<year>_<month>_(validated|interim)_reanalysis.nc
     ```

     con priorità `validated > interim`.

2. Per ogni combinazione **(anno, inquinante)** e per ogni regione in `REGIONS`:

   * normalizza le longitudini in `[-180, 180)`
   * effettua il clip sulla regione
   * salva la **mappa annuale** (PNG)
   * calcola e salva la **media regionale pesata** (CSV)
   * memorizza i frame per la GIF

3. Per ogni **(regione, inquinante)** crea:

   * **GIF** multi-anno
   * **montage multi-anno** (anni affiancati) con **colorbar unica**

4. *(Opzionale ma consigliato)* Estrae i valori nei capoluoghi e salva:

   * **time-series per capoluogo** (PNG)
   * **CSV dedicato ai capoluoghi**

---

## 3) Requisiti

* Python 3.x
* Pacchetti:

  * `numpy`, `pandas`, `xarray`
  * `matplotlib`
  * `geopandas`, `shapely`
  * `cartopy`
  * `imageio`
  * `rioxarray` *(opzionale ma consigliato)*
* Confini regionali:

  ```
  boundaries/regioni_italia.json
  ```
* Capoluoghi (una delle due opzioni):

  * `boundaries/capoluoghi_provincia.csv`
  * dizionario `PROVINCE_CAPITALS` nel codice

---

## 4) Config principali

Nel blocco **CONFIG**:

* `DATA_DIR` – cartella con i NetCDF CAMS
* `BOUNDARY_FILE` – GeoJSON regioni
* `OUT_DIR` – cartella output (`outputs_cams_annual`)
* `REGIONS` – regioni da processare (`Abruzzo`, `Puglia`)
* `PROVINCE_CAPITALS` – capoluoghi con lat/lon e sigla (AQ, TE, …)

---

## 5) Output: cosa viene creato e dove

### A) Mappe annuali (PNG)

**Percorso**

```text
outputs_cams_annual/<Regione>/<Anno>/
<Regione>_<Anno>_<pollutant>_annual_mean.png
```

**Contenuto**

* raster dell’inquinante
* contorno della regione
* griglia lat/lon con ticks
* **capoluoghi di provincia** (punto + sigla), se abilitati con:

  ```python
  plot_province_capitals(ax, region_name)
  ```

---

### B) CSV delle medie regionali pesate

**Percorso**

```text
outputs_cams_annual/regional_annual_means.csv
```

**Colonne**

* `region`
* `year`
* `pollutant`
* `value` (media regionale)
* `units`
* `var_name`
* `n_months`

La media è pesata con `cos(lat)` per approssimare l’area delle celle.

---

### C) GIF multi-anno

**Percorso**

```text
outputs_cams_annual/gifs/
<Regione>_<pollutant>_annual_mean.gif
```

Ogni frame è la mappa annuale già salvata.

---

### D) Montage multi-anno con colorbar unica

**Percorso**

```text
outputs_cams_annual/montages/<Regione>/
<Regione>_<pollutant>_<y0>-<y1>_montage.png
```

**Note importanti**

* `vmin` / `vmax` calcolati su **tutti gli anni** (percentili 2–98)
* colorbar su asse dedicato:

  * `fig.subplots_adjust(bottom=...)`
  * `cax = fig.add_axes([...])`
* salvataggio **senza** `bbox_inches="tight"`

---

### E) Time-series dei capoluoghi (PNG)

⚠️ **ATTENZIONE**
Le funzioni:

* `sample_at_point`
* `build_capitals_timeseries_rows`
* `plot_timeseries_by_pollutant`

**sono presenti**, ma **se non vengono chiamate in `main()` non viene salvato nulla**.

**Output consigliato**

```text
outputs_cams_annual/timeseries_capoluoghi/
<Regione>_<pollutant>_capitals_timeseries.png
```

**Logica di attivazione**

* durante il loop annuale, usa `da_year`
* accumula righe in `ts_rows`
* alla fine:

  * crea `df_ts`
  * salva CSV
  * salva i plot con:

    ```python
    plot_timeseries_by_pollutant(df_ts, TS_DIR)
    ```

---

## 6) Inquinanti trattati: descrizione, rischi e soglie

> **OMS (WHO)** = linee guida sanitarie
> **UE** = limiti normativi (Direttiva 2008/50/CE)

### PM2.5

* **Cos’è:** particolato fine ≤ 2.5 µm
* **Rischi:** cardiovascolari e respiratori, asma, BPCO

**OMS 2021**

* annua: **5 µg/m³**
* 24h: **15 µg/m³**

**UE**

* annua: **25 µg/m³**

---

### PM10

* **Cos’è:** particolato inalabile ≤ 10 µm
* **Rischi:** irritazione vie respiratorie

**OMS 2021**

* annua: **15 µg/m³**
* 24h: **45 µg/m³**

**UE**

* annua: **40 µg/m³**
* giornaliera: **50 µg/m³** (≤ 35 superamenti/anno)

---

### NO₂

* **Cos’è:** biossido di azoto (traffico, combustione)
* **Rischi:** asma, infezioni respiratorie

**OMS 2021**

* annua: **10 µg/m³**
* 24h: **25 µg/m³**

**UE**

* annua: **40 µg/m³**
* oraria: **200 µg/m³** (≤ 18 superamenti/anno)

---

### O₃

* **Cos’è:** ozono troposferico (secondario)
* **Rischi:** riduzione funzione polmonare

**OMS 2021**

* peak season: **60 µg/m³**

**UE**

* 8h: **120 µg/m³**
* ≤ 25 giorni/anno (media su 3 anni)

> Nota: per O₃ le **medie annue non sono la metrica più rappresentativa**.

---

## 7) Note pratiche (gotcha comuni)

### Mappe che “cambiano dimensione”

* causa: `bbox_inches="tight"`
* soluzione: **non usarlo** per le mappe singole

### Colorbar che invade il montage

* evitare `tight_layout()`
* usare asse dedicato `cax`
* salvare **senza** `bbox_inches="tight"`

### Capoluoghi non visibili

* funzione non chiamata
* `zorder` troppo basso
* font troppo piccolo
* padding troppo stretto (`pad_x`, `pad_y`)

---

## 8) Esecuzione

```bash
python cams_make_annual_maps.py
```

Verifica:

* mappe PNG
* CSV medie regionali
* GIF
* montage
* *(opzionale)* time-series capoluoghi

---

## 9) Disclaimer (importante)

* I prodotti CAMS sono ideali per **pattern spaziali e trend**, non sostituiscono le stazioni ufficiali.
* Le soglie OMS/UE vanno interpretate secondo la **metrica temporale corretta**.
* Le time-series per capoluogo rappresentano il **valore della griglia CAMS**, non una misura urbana reale “street-level”.

---
