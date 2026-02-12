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

In questo progetto usi i file NetCDF della **reanalisi regionale europea** di qualità dell’aria (mensili *validated* o *interim* e/o file annuali già pronti). Questi prodotti combinano modellistica e osservazioni e sono adatti per **analisi climatologiche, trend temporali e confronti spaziali**.

---

## 2) Problema noto: cambio di griglia dopo il 2019 (e perché succede)

In alcune serie CAMS può capitare che **la griglia spaziale cambi nel tempo** (risoluzione/allineamento/packaging del prodotto), ad esempio attorno al 2019–2020.

Effetti tipici:

* le mappe “cambiano pixel” tra un anno e l’altro
* le GIF sembrano “saltare”
* i campionamenti *nearest* in un punto (capoluoghi) cambiano cella non per variazione reale, ma per **cambio griglia**
* anche la media regionale può cambiare leggermente perché cambiano le celle che contribuiscono

✅ **Soluzione adottata nello script:** **regrid (interpolazione) su una griglia comune** *per ogni inquinante*, così tutti gli anni sono confrontabili sullo stesso supporto spaziale.

> Nota scientifica: per campi di **concentrazione media annua** e per analisi di **trend/pattern** il regrid su griglia comune è una procedura standard e normalmente preferibile rispetto a confrontare campi su griglie diverse.

---

## 3) Cosa fa lo script (in breve)

Lo script:

1. **Scansiona i file NetCDF** in `DATA_DIR`

   * Se trova file annuali:

     ```text
     annual_mean_<pollutant>_<year>.nc
     ```

     li usa direttamente (più veloce).

   * Altrimenti calcola la media annua dai mensili:

     ```text
     cams-europe-air-quality-reanalyses_<pollutant>_<year>_<month>_(validated|interim)_reanalysis.nc
     ```

     con priorità `validated > interim`.

2. Per ogni combinazione **(anno, inquinante)**:

   * normalizza le longitudini in `[-180, 180)`
   * **imposta una griglia di riferimento** per quell’inquinante (vedi sezione “Regrid”)
   * **regrida tutti gli anni sulla griglia di riferimento**
   * calcola le **time-series dei capoluoghi** su griglia regridded (coerente nel tempo)
   * fa il **clip** per regione
   * salva la **mappa annuale** (PNG)
   * calcola e salva la **media regionale pesata** (CSV)
   * memorizza i frame per la GIF

3. Per ogni **(regione, inquinante)** crea:

   * **GIF** multi-anno
   * **montage multi-anno** (anni affiancati) con **colorbar unica**
     (anche qui i dati sono già regriddati prima del clip)

---

## 4) Regrid: cosa fa e come si configura

### Perché si fa

Il regrid serve a garantire che:

* ogni anno abbia **la stessa griglia**
* GIF e montage siano coerenti
* i campionamenti nei capoluoghi non “saltino” per cambio cella
* i confronti interannuali siano più solidi

### Come funziona

* si definisce una **griglia di riferimento** per ogni `pollutant`
* ogni `da_year` viene interpolato su quella griglia tramite `xarray.interp`

### Parametri principali (CONFIG)

* `REGRID_MODE`

  * `"first"`: usa la **prima griglia incontrata** per quel pollutant (consigliato: più stabile e riproducibile)
  * `"finest"`: usa la **griglia più fine** tra quelle incontrate (più dettaglio, più costo)
* `REGRID_METHOD`

  * `"linear"`: bilineare (default, consigliato per mappe/trend)
  * `"nearest"`: più “a blocchi”, ma evita smoothing

---

## 5) Requisiti

* Python 3.x
* Pacchetti:

  * `numpy`, `pandas`, `xarray`
  * `matplotlib`
  * `geopandas`, `shapely`
  * `cartopy`
  * `imageio`
  * `rioxarray` *(opzionale ma consigliato)*

File necessari:

* Confini regionali:

  ```text
  boundaries/regioni_italia.json
  ```
* Capoluoghi:

  * `boundaries/capoluoghi_provincia.csv` *(opzionale; nel codice c’è già `PROVINCE_CAPITALS`)*

---

## 6) Config principali

Nel blocco **CONFIG**:

* `DATA_DIR` – cartella con i NetCDF CAMS
* `BOUNDARY_FILE` – GeoJSON regioni
* `OUT_DIR` – cartella output (`outputs_cams_annual`)
* `REGIONS` – regioni da processare (`Abruzzo`, `Puglia`)
* `PROVINCE_CAPITALS` – capoluoghi con lat/lon e sigla (AQ, TE, …)
* `REGRID_MODE`, `REGRID_METHOD` – controlli per griglia comune

---

## 7) Output: cosa viene creato e dove

### A) Mappe annuali (PNG)

Percorso:

```text
outputs_cams_annual/<Regione>/<Anno>/
<Regione>_<Anno>_<pollutant>_annual_mean.png
```

Contenuto:

* raster dell’inquinante (**regriddato** su griglia comune)
* contorno della regione
* griglia lat/lon con ticks
* capoluoghi di provincia (punto + sigla)

---

### B) CSV delle medie regionali pesate

Percorso:

```text
outputs_cams_annual/regional_annual_means.csv
```

Colonne:

* `region`
* `year`
* `pollutant`
* `value` (media regionale)
* `units`
* `var_name`
* `n_months`

La media è pesata con `cos(lat)` per approssimare l’area delle celle.

> Nota: essendo su griglia comune, la comparabilità temporale della media regionale è più stabile.

---

### C) GIF multi-anno

Percorso:

```text
outputs_cams_annual/gifs/
<Regione>_<pollutant>_annual_mean.gif
```

Ogni frame è la mappa annuale già salvata.

---

### D) Montage multi-anno con colorbar unica

Percorso:

```text
outputs_cams_annual/montages/<Regione>/
<Regione>_<pollutant>_<y0>-<y1>_montage.png
```

Note:

* `vmin/vmax` calcolati su tutti gli anni (percentili 2–98)
* colorbar su asse dedicato (`cax = fig.add_axes([...])`)
* salvataggio senza `bbox_inches="tight"`

---

### E) Time-series dei capoluoghi (PNG + CSV)

Percorso:

```text
outputs_cams_annual/timeseries_capoluoghi/
<Regione>_<pollutant>_capitals_timeseries.png
outputs_cams_annual/timeseries_capoluoghi/capitals_timeseries_values.csv
```

Importante:

* i valori sono estratti su **griglia regridded** → evita salti artificiali dovuti a cambio griglia
* rappresentano il valore della cella CAMS, non una misura “street-level”

---

## 8) Inquinanti trattati: descrizione, rischi e soglie

> OMS (WHO) = linee guida sanitarie
> UE = limiti normativi (Direttiva 2008/50/CE)

### PM2.5

* OMS 2021: annua **5 µg/m³**, 24h **15 µg/m³**
* UE: annua **25 µg/m³**

### PM10

* OMS 2021: annua **15 µg/m³**, 24h **45 µg/m³**
* UE: annua **40 µg/m³**, 24h **50 µg/m³** (≤ 35 superamenti/anno)

### NO₂

* OMS 2021: annua **10 µg/m³**, 24h **25 µg/m³**
* UE: annua **40 µg/m³**, oraria **200 µg/m³** (≤ 18 superamenti/anno)

### O₃

* OMS 2021: peak season **60 µg/m³**
* UE: 8h **120 µg/m³** (≤ 25 giorni/anno, media su 3 anni)

> Nota: per O₃, la media annua non è la metrica più rappresentativa.

---

## 9) Note pratiche (gotcha comuni)

### Mappe che “cambiano” nel tempo

* causa: cambio griglia nel dataset (o mix annual/mensili, validated/interim)
* soluzione: **regrid su griglia comune** (attivo nello script)

### Mappe che cambiano dimensione/impaginazione

* causa: `bbox_inches="tight"`
* soluzione: non usarlo sulle mappe singole (già fatto)

### Capoluoghi non visibili

* controllare che `plot_province_capitals(...)` sia chiamata
* `zorder`/font/padding

---

## 10) Esecuzione

```bash
python cams_make_annual_maps.py
```

Verifica:

* mappe PNG
* CSV medie regionali
* GIF
* montage
* time-series capoluoghi

---

## 11) Disclaimer (importante)

* I prodotti CAMS sono ideali per **pattern spaziali e trend**, non sostituiscono le stazioni ufficiali.
* Le soglie OMS/UE vanno interpretate con la **metrica temporale corretta**.
* Le time-series per capoluogo rappresentano il **valore della griglia CAMS**.
* Il **regrid** è usato per rendere i confronti interannuali coerenti; va dichiarato se usi output in report/paper.


