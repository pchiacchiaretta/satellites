# CAMS – Mappe annuali, medie regionali, GIF, montage e time-series per capoluoghi (Abruzzo + Puglia)

Questo progetto genera **mappe annuali** di inquinanti atmosferici (da CAMS), calcola **medie regionali**, crea **GIF** nel tempo, costruisce un **montage multi-anno** (anni affiancati con **un’unica colorbar**) e può produrre **grafici di andamento temporale (time-series) per capoluogo di provincia**.

---

## 1) Cos’è CAMS e cosa stai usando qui

**CAMS (Copernicus Atmosphere Monitoring Service)** è il servizio Copernicus dedicato al monitoraggio atmosferico: produce e distribuisce **analisi**, **previsioni** e **reanalisi** (ricostruzioni “storiche” consistenti) per vari componenti atmosferici e qualità dell’aria. :contentReference[oaicite:0]{index=0}

Nel tuo caso usi i file NetCDF della **reanalisi/produzione regionale europea** di qualità dell’aria (dataset mensili “validated” o “interim” e/o file già annuali). Questi prodotti tipicamente combinano modellistica e osservazioni, e sono pensati per analisi climatologiche/temporali e confronti tra aree. :contentReference[oaicite:1]{index=1}

---

## 2) Cosa fa lo script (in breve)

Lo script:

1. **Scansiona i file NetCDF** in `DATA_DIR`  
   - Se trova file annuali: `annual_mean_<pollutant>_<year>.nc` li usa direttamente (più veloce).
   - Altrimenti calcola la media annua dai mensili:
     `cams-europe-air-quality-reanalyses_<pollutant>_<year>_<month>_(validated|interim)_reanalysis.nc`
     con priorità `validated > interim`.

2. Per ogni combinazione (anno, inquinante) e per ogni regione in `REGIONS`:
   - normalizza le longitudini in `[-180,180)`,
   - fa clip sulla regione,
   - salva **mappa annuale** (PNG),
   - calcola e salva **media regionale pesata** (CSV),
   - memorizza i frame per costruire la GIF.

3. Per ogni (regione, inquinante) crea:
   - **GIF** multi-anno.
   - **Montage** multi-anno (anni affiancati) con **colorbar unica**.

4. (Opzionale, ma consigliato) Estrae valori ai capoluoghi e salva:
   - **time-series per capoluogo** (un PNG per ogni coppia regione-inquinante),
   - e/o un CSV dedicato ai capoluoghi.

---

## 3) Requisiti

- Python 3.x
- pacchetti: `numpy`, `pandas`, `xarray`, `matplotlib`, `geopandas`, `shapely`, `cartopy`, `imageio`, `rioxarray` (opzionale ma consigliato)
- file di confine regioni: `boundaries/regioni_italia.json`
- (opzionale) file capoluoghi: `boundaries/capoluoghi_provincia.csv` **oppure** dizionario `PROVINCE_CAPITALS` già nel codice

---

## 4) Config principali

Nel blocco CONFIG:

- `DATA_DIR`: cartella con i NetCDF CAMS
- `BOUNDARY_FILE`: GeoJSON regioni
- `OUT_DIR`: cartella output (`outputs_cams_annual`)
- `REGIONS`: elenco regioni da processare (Abruzzo, Puglia)
- `PROVINCE_CAPITALS`: elenco capoluoghi con lat/lon + sigla (AQ, TE, …)

---

## 5) Output: cosa viene creato e dove

### A) Mappe annuali per regione/anno/inquinante (PNG)

Percorso:
```

outputs_cams_annual/<Regione>/<Anno>/<Regione>*<Anno>*<pollutant>_annual_mean.png

````

Contenuto:
- raster dell’inquinante
- contorno regione
- griglia lat/lon e ticks
- **capoluoghi di provincia** (punto + sigla) se abilitati in `plot_map()` con:
  ```python
  plot_province_capitals(ax, region_name)
````

---

### B) CSV delle medie regionali pesate

Percorso:

```
outputs_cams_annual/regional_annual_means.csv
```

Colonne tipiche:

* `region`, `year`, `pollutant`
* `value` (media regionale)
* `units`, `var_name`, `n_months`

La media è pesata con `cos(lat)` per approssimare l’area delle celle.

---

### C) GIF (una per regione+inquinante)

Percorso:

```
outputs_cams_annual/gifs/<Regione>_<pollutant>_annual_mean.gif
```

Ogni frame è la mappa annuale già salvata.

---

### D) Montage multi-anno (anni affiancati) con colorbar unica

Percorso:

```
outputs_cams_annual/montages/<Regione>/<Regione>_<pollutant>_<y0>-<y1>_montage.png
```

Note importanti:

* vmin/vmax calcolati su tutti gli anni (percentili 2–98) per avere scala coerente.
* Per evitare che la colorbar “invada” i subplot, la funzione migliore è quella dove:

  * fai `fig.subplots_adjust(..., bottom=...)`
  * crei un asse dedicato `cax = fig.add_axes([...])`
  * salvi **senza** `bbox_inches="tight"`.

---

### E) Time-series capoluoghi (PNG) – **ATTENZIONE: va chiamata in `main()`**

Le funzioni `sample_at_point`, `build_capitals_timeseries_rows`, `plot_timeseries_by_pollutant`
nel tuo codice **ci sono**, ma **se non le chiami in `main()` non verrà salvato nulla**.

Consiglio output:

Percorso:

```
outputs_cams_annual/timeseries_capoluoghi/<Regione>_<pollutant>_capitals_timeseries.png
```

**Dove lo salva?**
Se imposti, ad esempio:

```python
TS_DIR = OUT_DIR / "timeseries_capoluoghi"
plot_timeseries_by_pollutant(df_ts, TS_DIR)
```

allora i grafici finiscono lì.

✅ **Come abilitarlo (logica)**
Durante il loop principale, dopo aver creato `da_year` (annuale) e prima/dopo il loop regioni, costruisci righe time-series per capoluoghi (usando `da_year` già normalizzato) e accumuli in una lista `ts_rows`. Alla fine:

* `df_ts = pd.DataFrame(ts_rows)`
* salvi CSV
* salvi i plot con `plot_timeseries_by_pollutant(df_ts, TS_DIR)`

---

## 6) Inquinanti (“composti”) trattati: cosa sono, rischi salute, soglie

Qui sotto trovi una sintesi **sanitaria** e i riferimenti a **valori guida OMS 2021** e **limiti UE (Direttiva 2008/50/CE)**.

> Nota: OMS (WHO) = **linee guida** sanitarie (più stringenti); UE = **limiti normativi** (obblighi legali). Le metriche possono differire (annuale, giornaliero, 8h, ecc.).

### PM2.5 (particulate_matter_2.5um)

**Cos’è:** particolato fine con diametro ≤ 2.5 µm.
**Perché è dannoso:** penetra in profondità nei polmoni e può entrare nel circolo sanguigno; associato a mortalità cardiovascolare/respiratoria e peggioramento asma/BPCO.

**Valori guida OMS 2021 (PM2.5):**

* Media annua: **5 µg/m³**
* Media 24h: **15 µg/m³**

**Limiti UE (Direttiva 2008/50/CE):**

* Media annua: **25 µg/m³** (valore limite) ([confluence.ecmwf.int][1])

---

### PM10 (particulate_matter_10um)

**Cos’è:** particolato inalabile con diametro ≤ 10 µm.
**Perché è dannoso:** irritazione vie respiratorie, peggioramento patologie respiratorie; parte può depositarsi in bronchi e polmoni.

**Valori guida OMS 2021 (PM10):**

* Media annua: **15 µg/m³**
* Media 24h: **45 µg/m³**

**Limiti UE (Direttiva 2008/50/CE):**

* Media annua: **40 µg/m³**
* Giornaliero: **50 µg/m³** da non superare > **35** volte/anno ([confluence.ecmwf.int][1])

---

### NO2 (nitrogen_dioxide)

**Cos’è:** biossido di azoto, prodotto soprattutto da combustione (traffico, riscaldamento).
**Perché è dannoso:** irritazione delle vie aeree, aumento suscettibilità a infezioni, peggioramento asma; è anche indicatore di inquinamento da traffico.

**Valori guida OMS 2021 (NO2):**

* Media annua: **10 µg/m³**
* Media 24h: **25 µg/m³**

**Limiti UE (Direttiva 2008/50/CE):**

* Media annua: **40 µg/m³**
* Orario: **200 µg/m³** da non superare > **18** volte/anno ([confluence.ecmwf.int][1])

---

### O3 (ozone)

**Cos’è:** ozono troposferico (inquinante “secondario”), si forma con reazioni fotochimiche da precursori (NOx, VOC) con sole.
**Perché è dannoso:** irritazione, riduzione funzione polmonare, peggioramento asma; effetti più marcati nei picchi estivi.

**Valori guida OMS 2021 (O3):**

* Picco stagione: **60 µg/m³** (metrica OMS “peak season”)

**Soglie UE (Direttiva 2008/50/CE):**

* Valore obiettivo per la protezione della salute: **120 µg/m³** (massima media mobile 8 ore), da non superare > **25** giorni/anno (media su 3 anni) ([confluence.ecmwf.int][1])

> Importante: per O3 le “medie annue” non sono la metrica più rappresentativa dei rischi: contano molto i picchi (8h, stagionali).

---

## 7) Note pratiche e “gotcha” comuni

### Perché a volte la mappa “cambiava dimensione”?

Di solito succede quando si salva con `bbox_inches="tight"`: Matplotlib ritaglia diversamente in base a testi/ticks, e la GIF “balla”.
Nel tuo `plot_map()` fai bene: **non usare** `bbox_inches="tight"` per le mappe singole.

### Perché nel montage la colorbar finiva sopra i subplot?

Perché `tight_layout()` + `bbox_inches="tight"` possono ricalcolare i bounding box e “tirare su” la colorbar.
Soluzione corretta (che hai già iniziato a usare): **asse dedicato `cax`** + `subplots_adjust(bottom=...)` + salvataggio **senza** `bbox_inches="tight"`.

### Capoluoghi non visibili sulle mappe

Cause tipiche:

* chiamata a `plot_province_capitals(ax, region_name)` mancante (ora l’hai messa)
* zorder troppo basso (tu usi 10/11, va bene)
* label troppo piccola o coperta: aumenta `fontsize` o `bbox alpha`
* coordinate leggermente fuori extent per padding troppo stretto: aumenta `pad_x/pad_y` (es. 0.05 invece di 0.03)

---

## 8) Esecuzione

Esegui:

```bash
python cams_make_annual_maps.py
```

Controlla poi:

* `outputs_cams_annual/<Regione>/<Anno>/...png` (mappe)
* `outputs_cams_annual/regional_annual_means.csv`
* `outputs_cams_annual/gifs/...gif`
* `outputs_cams_annual/montages/...png`
* (se abilitato) `outputs_cams_annual/timeseries_capoluoghi/...png`

---

## 9) Disclaimer (molto importante)

* Le reanalisi/modelli CAMS sono ottimi per **pattern spaziali e trend**, ma non sostituiscono misure ufficiali puntuali di stazioni locali.
* Le soglie OMS/UE vanno interpretate correttamente rispetto alla **metrica temporale** (annuale, 24h, 8h, oraria).
* I grafici “per capoluogo” estraggono il valore della **griglia CAMS nel punto** (nearest/interp), non una misura reale urbana “street level”.

