#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CAMS – Mappe annuali + CSV medie regionali + GIF per (regione, inquinante) nel tempo.

- Legge:
  1) annual_mean_<pollutant>_<year>.nc (se presenti)  [più veloce e stabile]
  2) altrimenti mensili cams-europe-air-quality-reanalyses_<pollutant>_<year>_<month>_(validated|interim)_reanalysis.nc

- Priorità mensili: validated > interim

FIX (importante):
- Regrid su griglia comune (per ogni pollutant) per evitare cambi di griglia dopo 2019 e rendere confronti/GIF coerenti.
"""

import os
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"  # filesystem condivisi
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

import geopandas as gpd
from shapely.geometry import mapping

import imageio.v2 as imageio

import cartopy.crs as ccrs
import cartopy.feature as cfeature

import matplotlib.ticker as mticker
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter


# =========================
# CONFIG
# =========================
DATA_DIR = Path("/home/manager/datastorage/atmolab/data/projects/onco_mammo/satellite/cams")
BOUNDARY_FILE = Path("boundaries/regioni_italia.json")

OUT_DIR = Path("outputs_cams_annual")
GIF_DIR = OUT_DIR / "gifs"

REGIONS = ["Abruzzo", "Puglia"]

DEFAULT_CMAP = "viridis"
DPI = 200

FOOTER_LEFT = "Dati: Copernicus Atmosphere Monitoring Service (CAMS) – Copernicus/ECMWF"
FOOTER_RIGHT = "Elaborazione: UdAtmo - Dr. Piero Chiacchiaretta"

GIF_FPS = 1.2
GIF_LOOP = 0

CAPOLUOGHI_FILE = Path("boundaries/capoluoghi_provincia.csv")

# Regrid config:
# - "first": usa la prima griglia incontrata per quel pollutant (più riproducibile)
# - "finest": usa la griglia con passo più fine tra quelle incontrate (più dettagli, più costoso)
REGRID_MODE = "first"   # ["first", "finest"]
REGRID_METHOD = "linear"  # ["linear", "nearest"]
DEBUG_GRID = False        # True per print su shape/step


df_caps = pd.read_csv(CAPOLUOGHI_FILE)
gdf_caps = gpd.GeoDataFrame(
    df_caps,
    geometry=gpd.points_from_xy(df_caps.lon, df_caps.lat),
    crs="EPSG:4326"
)

PROVINCE_CAPITALS = {
    "Abruzzo": [
        {"name": "L'Aquila", "sigla": "AQ", "lat": 42.3496, "lon": 13.3995},
        {"name": "Teramo",   "sigla": "TE", "lat": 42.6589, "lon": 13.7044},
        {"name": "Pescara",  "sigla": "PE", "lat": 42.4618, "lon": 14.2161},
        {"name": "Chieti",   "sigla": "CH", "lat": 42.3510, "lon": 14.1675},
    ],
    "Puglia": [
        {"name": "Bari",      "sigla": "BA", "lat": 41.1171, "lon": 16.8719},
        {"name": "Brindisi",  "sigla": "BR", "lat": 40.6380, "lon": 17.9450},
        {"name": "Foggia",    "sigla": "FG", "lat": 41.4622, "lon": 15.5446},
        {"name": "Lecce",     "sigla": "LE", "lat": 40.3529, "lon": 18.1743},
        {"name": "Taranto",   "sigla": "TA", "lat": 40.4644, "lon": 17.2470},
        {"name": "Barletta",  "sigla": "BT", "lat": 41.3207, "lon": 16.2810},
    ],
}


# =========================
# REGEX FILE
# =========================
ANNUAL_RE = re.compile(
    r"annual_mean_(?P<pollutant>.+?)_(?P<year>\d{4})\.nc$",
    re.IGNORECASE
)

MONTHLY_RE = re.compile(
    r"cams-europe-air-quality-reanalyses_(?P<pollutant>.+?)_(?P<year>\d{4})_(?P<month>\d{2})_(?P<kind>interim_reanalysis|validated_reanalysis)\.nc$",
    re.IGNORECASE
)


# =========================
# UTILS: lat/lon + main var
# =========================
def find_lat_lon(ds: xr.Dataset) -> Tuple[str, str]:
    lat_candidates = ["latitude", "lat", "Latitude", "nav_lat"]
    lon_candidates = ["longitude", "lon", "Longitude", "nav_lon"]

    lat = next((c for c in lat_candidates if c in ds.coords), None)
    lon = next((c for c in lon_candidates if c in ds.coords), None)

    if lat is None:
        lat = next((v for v in lat_candidates if v in ds.variables), None)
    if lon is None:
        lon = next((v for v in lon_candidates if v in ds.variables), None)

    if lat is None or lon is None:
        raise ValueError("Impossibile trovare coordinate lat/lon nel dataset.")
    return lat, lon


def guess_main_data_var(ds: xr.Dataset) -> str:
    blacklist = set(ds.coords.keys()) | {"crs"}
    candidates = []
    for v in ds.data_vars:
        if v in blacklist:
            continue
        da = ds[v]
        candidates.append((v, da.ndim, da.size))
    if not candidates:
        raise ValueError("Nessuna data_var trovata per il plot.")
    candidates.sort(key=lambda x: (x[1], x[2]), reverse=True)
    return candidates[0][0]


def normalize_longitudes(da: xr.DataArray, lon_name: str) -> xr.DataArray:
    lon = da[lon_name]
    lon_norm = ((lon + 180) % 360) - 180
    da2 = da.assign_coords({lon_name: lon_norm})
    return da2.sortby(lon_name)


def _grid_step_1d(coord: np.ndarray) -> float:
    coord = np.asarray(coord)
    if coord.size < 2:
        return np.nan
    d = np.diff(coord.astype(float))
    return float(np.nanmedian(np.abs(d)))


# =========================
# REGRID (NUOVO)
# =========================
def regrid_to_reference(
    da: xr.DataArray,
    ref_lats: xr.DataArray,
    ref_lons: xr.DataArray,
    method: str = "linear",
) -> xr.DataArray:
    """
    Regrid su griglia (ref_lats, ref_lons) usando xarray.interp.
    """
    lat_name = da.attrs["__lat_name__"]
    lon_name = da.attrs["__lon_name__"]

    # lon coerenti e ordinate
    da = normalize_longitudes(da, lon_name)

    da_rg = da.interp(
        {lat_name: ref_lats, lon_name: ref_lons},
        method=method,
        kwargs={"fill_value": np.nan},
    )
    da_rg.attrs.update(da.attrs)
    return da_rg


def maybe_set_reference_grid(
    pollutant: str,
    da: xr.DataArray,
    ref_grid_by_pollutant: Dict[str, Tuple[xr.DataArray, xr.DataArray]],
    mode: str = "first",
) -> None:
    """
    Imposta (se necessario) la griglia di riferimento per quel pollutant.
    mode:
      - "first": prima griglia incontrata
      - "finest": prende la più fine (min step) tra quelle incontrate
    """
    lat_name = da.attrs["__lat_name__"]
    lon_name = da.attrs["__lon_name__"]

    lats = da[lat_name]
    lons = da[lon_name]

    if pollutant not in ref_grid_by_pollutant:
        ref_grid_by_pollutant[pollutant] = (lats.copy(), lons.copy())
        return

    if mode != "finest":
        return

    ref_lats, ref_lons = ref_grid_by_pollutant[pollutant]
    step_lat_new = _grid_step_1d(lats.values)
    step_lon_new = _grid_step_1d(lons.values)
    step_lat_ref = _grid_step_1d(ref_lats.values)
    step_lon_ref = _grid_step_1d(ref_lons.values)

    # "finest" = step più piccolo (più risoluzione)
    if (np.isfinite(step_lat_new) and np.isfinite(step_lat_ref) and step_lat_new < step_lat_ref) or \
       (np.isfinite(step_lon_new) and np.isfinite(step_lon_ref) and step_lon_new < step_lon_ref):
        ref_grid_by_pollutant[pollutant] = (lats.copy(), lons.copy())


# =========================
# CONFINE REGIONI
# =========================
def load_region_geoms() -> Dict[str, object]:
    gdf = gpd.read_file(BOUNDARY_FILE)
    region_geoms = {}
    for r in REGIONS:
        mask = gdf["name"].astype(str).str.lower().eq(r.lower())
        reg = gdf.loc[mask]
        if reg.empty:
            raise ValueError(f"Regione '{r}' non trovata nel GeoJSON.")
        region_geoms[r] = reg.geometry.iloc[0]
    return region_geoms


# =========================
# SCANSIONE FILE
# =========================
def scan_cams_files(data_dir: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    annual_rows = []
    monthly_rows = []

    for f in data_dir.rglob("*.nc"):
        name = f.name

        ma = ANNUAL_RE.search(name)
        if ma:
            annual_rows.append((
                int(ma.group("year")),
                ma.group("pollutant").lower(),
                str(f)
            ))
            continue

        mm = MONTHLY_RE.search(name)
        if mm:
            monthly_rows.append((
                int(mm.group("year")),
                int(mm.group("month")),
                mm.group("pollutant").lower(),
                mm.group("kind").lower(),
                str(f)
            ))

    df_annual = pd.DataFrame(annual_rows, columns=["year", "pollutant", "path"]).sort_values(["year", "pollutant"])
    df_monthly = pd.DataFrame(monthly_rows, columns=["year", "month", "pollutant", "kind", "path"]).sort_values(
        ["year", "pollutant", "month", "kind"]
    )

    return df_annual, df_monthly


# =========================
# CALCOLO MEDIA ANNUA
# =========================
def compute_annual_from_single_nc(path: str) -> xr.DataArray:
    with xr.open_dataset(path) as ds:
        lat_name, lon_name = find_lat_lon(ds)
        var_name = guess_main_data_var(ds)

        da = ds[var_name]

        time_dim = next((d for d in da.dims if d.lower() in ("time", "valid_time", "datetime")), None)
        if time_dim is not None:
            da = da.mean(time_dim, skipna=True)

        da = da.astype("float64")

        attrs = dict(da.attrs)
        if not attrs.get("units"):
            attrs["units"] = ds[var_name].attrs.get("units", "") or ds.attrs.get("units", "")
        if not attrs.get("long_name"):
            attrs["long_name"] = ds[var_name].attrs.get("long_name", "") or ds.attrs.get("title", "")

        da.attrs.update(attrs)
        da.attrs["__lat_name__"] = lat_name
        da.attrs["__lon_name__"] = lon_name
        da.attrs["__var_name__"] = var_name
        da.attrs["__n_months__"] = 12
        return da


def compute_yearly_mean_from_monthlies(paths: List[str]) -> xr.DataArray:
    paths = sorted(paths)

    da_sum = None
    n_months = 0

    lat_name = lon_name = var_name = None

    for p in paths:
        with xr.open_dataset(p) as ds:
            if lat_name is None or lon_name is None:
                lat_name, lon_name = find_lat_lon(ds)
            if var_name is None:
                var_name = guess_main_data_var(ds)

            da = ds[var_name]

            time_dim = next((d for d in da.dims if d.lower() in ("time", "valid_time", "datetime")), None)
            if time_dim is not None:
                da_m = da.mean(time_dim, skipna=True)
            else:
                da_m = da.mean(da.dims[0], skipna=True) if da.ndim > 2 else da

            da_m = da_m.astype("float64")

            if da_sum is None:
                da_sum = da_m
                da_sum.attrs.update(da.attrs)

                if not da_sum.attrs.get("units"):
                    da_sum.attrs["units"] = ds[var_name].attrs.get("units", "") or ds.attrs.get("units", "")
                if not da_sum.attrs.get("long_name"):
                    da_sum.attrs["long_name"] = ds[var_name].attrs.get("long_name", "") or ds.attrs.get("title", "")
            else:
                da_sum = da_sum + da_m

            n_months += 1

    if da_sum is None or n_months == 0:
        raise RuntimeError("Nessun file mensile valido per calcolare la media annua.")

    da_y = da_sum / n_months
    da_y.attrs.update(da_sum.attrs)
    da_y.attrs["__lat_name__"] = lat_name
    da_y.attrs["__lon_name__"] = lon_name
    da_y.attrs["__var_name__"] = var_name
    da_y.attrs["__n_months__"] = n_months
    return da_y


def pick_monthly_paths_for_year_pollutant(grp: pd.DataFrame) -> List[str]:
    grp_v = grp[grp["kind"].str.contains("validated")]
    grp_i = grp[grp["kind"].str.contains("interim")]
    use = grp_v if len(grp_v) > 0 else grp_i
    return use["path"].tolist()


# =========================
# CLIP + MEDIA REGIONALE
# =========================
def clip_to_region(da: xr.DataArray, lat_name: str, lon_name: str, geom):
    minx, miny, maxx, maxy = geom.bounds

    da_sub = da.where(
        (da[lat_name] >= miny) & (da[lat_name] <= maxy) &
        (da[lon_name] >= minx) & (da[lon_name] <= maxx),
        drop=False
    )

    try:
        import rioxarray  # noqa: F401
        da_sub = da_sub.rio.write_crs("EPSG:4326", inplace=False)
        da_clip = da_sub.rio.clip([mapping(geom)], crs="EPSG:4326", drop=False)
        return da_clip
    except Exception:
        return da_sub


def regional_weighted_mean(da_reg: xr.DataArray) -> float:
    lat_name = da_reg.attrs["__lat_name__"]
    lat_radians = np.deg2rad(da_reg[lat_name])
    weights = np.cos(lat_radians)
    wmean = da_reg.weighted(weights).mean(skipna=True)
    return float(wmean.values)


# =========================
# PLOT
# =========================
def plot_province_capitals(ax, region_name: str):
    if region_name not in PROVINCE_CAPITALS:
        return

    for cap in PROVINCE_CAPITALS[region_name]:
        lon = cap["lon"]
        lat = cap["lat"]
        sigla = cap["sigla"]

        ax.plot(
            lon, lat,
            marker="o",
            markersize=4,
            color="black",
            transform=ccrs.PlateCarree(),
            zorder=10,
        )

        ax.text(
            lon + 0.03, lat + 0.03,
            sigla,
            fontsize=9,
            fontweight="bold",
            ha="left",
            va="bottom",
            transform=ccrs.PlateCarree(),
            zorder=11,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.6, pad=0.8),
        )


def plot_map(
    da_reg: xr.DataArray,
    region_name: str,
    year: int,
    pollutant: str,
    geom,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
) -> Path:
    lat_name = da_reg.attrs["__lat_name__"]
    lon_name = da_reg.attrs["__lon_name__"]

    units = da_reg.attrs.get("units") or "µg/m3"
    long_name = (da_reg.attrs.get("long_name") or da_reg.attrs.get("standard_name") or pollutant)

    out_dir = OUT_DIR / region_name / str(year)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / f"{region_name}_{year}_{pollutant}_annual_mean.png"

    lats = np.asarray(da_reg[lat_name].values)
    lons = np.asarray(da_reg[lon_name].values)
    data = np.ma.masked_invalid(np.asarray(da_reg.values))

    if data.size == 0 or lats.size == 0 or lons.size == 0:
        raise ValueError("Griglia vuota dopo clip (lat/lon/data).")
    if np.ma.count(data) == 0:
        raise ValueError("Tutti i valori sono NaN dopo clip.")

    minx, miny, maxx, maxy = geom.bounds
    pad_x = (maxx - minx) * 0.03
    pad_y = (maxy - miny) * 0.03
    extent = [minx - pad_x, maxx + pad_x, miny - pad_y, maxy + pad_y]

    fig = plt.figure(figsize=(8.6, 7.6), dpi=DPI)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(extent, crs=ccrs.PlateCarree())

    xmaj = np.arange(np.floor(extent[0]), np.ceil(extent[1]) + 1e-6, 1.0)
    ymaj = np.arange(np.floor(extent[2]), np.ceil(extent[3]) + 1e-6, 0.5)
    xmin = np.arange(np.floor(extent[0]), np.ceil(extent[1]) + 1e-6, 0.5)
    ymin = np.arange(np.floor(extent[2]), np.ceil(extent[3]) + 1e-6, 0.25)

    ax.set_xticks(xmaj, crs=ccrs.PlateCarree())
    ax.set_yticks(ymaj, crs=ccrs.PlateCarree())
    ax.xaxis.set_minor_locator(mticker.FixedLocator(xmin))
    ax.yaxis.set_minor_locator(mticker.FixedLocator(ymin))

    ax.xaxis.set_major_formatter(LongitudeFormatter(number_format=".0f", degree_symbol="°"))
    ax.yaxis.set_major_formatter(LatitudeFormatter(number_format=".1f", degree_symbol="°"))

    ax.tick_params(axis="both", which="major", direction="in", top=True, right=True, length=7, labelsize=10)
    ax.tick_params(axis="both", which="minor", direction="in", top=True, right=True, length=4)
    ax.tick_params(axis="x", labeltop=True, labelbottom=True)
    ax.tick_params(axis="y", labelright=True, labelleft=True)

    ax.grid(True, which="major", linewidth=0.4, linestyle="--", alpha=0.5)
    ax.grid(True, which="minor", linewidth=0.25, linestyle=":", alpha=0.35)

    ax.add_feature(cfeature.COASTLINE.with_scale("10m"), linewidth=0.6, zorder=2)
    ax.add_feature(cfeature.BORDERS.with_scale("10m"), linewidth=0.4, zorder=2)

    mesh = ax.pcolormesh(
        lons, lats, data,
        transform=ccrs.PlateCarree(),
        shading="auto",
        cmap=DEFAULT_CMAP,
        vmin=vmin,
        vmax=vmax,
        zorder=1,
    )

    ax.add_geometries([geom], crs=ccrs.PlateCarree(),
                      facecolor="none", edgecolor="black", linewidth=1.8, zorder=5)

    plot_province_capitals(ax, region_name)
    ax.set_title(f"{region_name} – Media annua {year}\n{pollutant}", fontsize=13)

    cb = plt.colorbar(mesh, ax=ax, orientation="horizontal", pad=0.08, fraction=0.06, aspect=35)
    cb.set_label(f"{long_name} [{units}]" if units else f"{long_name}", fontsize=10)
    cb.ax.tick_params(labelsize=9, direction="in")

    fig.text(0.01, 0.01, FOOTER_LEFT, ha="left", va="bottom", fontsize=7)
    fig.text(0.99, 0.01, FOOTER_RIGHT, ha="right", va="bottom", fontsize=7)

    fig.savefig(out_png, dpi=DPI)
    plt.close(fig)
    return out_png


# =========================
# TIME SERIES (capoluoghi)
# =========================
def sample_at_point(da: xr.DataArray, lat: float, lon: float, method: str = "nearest") -> float:
    lat_name = da.attrs["__lat_name__"]
    lon_name = da.attrs["__lon_name__"]
    lon_norm = ((lon + 180) % 360) - 180

    if method == "nearest":
        v = da.sel({lat_name: lat, lon_name: lon_norm}, method="nearest")
    else:
        v = da.interp({lat_name: lat, lon_name: lon_norm}, method="linear")

    return float(v.values)


def build_capitals_timeseries_rows(
    da_year: xr.DataArray,
    year: int,
    pollutant: str,
    region_name: str,
    method: str = "nearest",
):
    units = da_year.attrs.get("units", "") or "µg/m3"
    long_name = da_year.attrs.get("long_name") or da_year.attrs.get("standard_name") or pollutant

    rows = []
    for cap in PROVINCE_CAPITALS.get(region_name, []):
        try:
            val = sample_at_point(da_year, cap["lat"], cap["lon"], method=method)
        except Exception:
            val = np.nan

        rows.append({
            "region": region_name,
            "capital": cap["name"],
            "sigla": cap["sigla"],
            "year": year,
            "pollutant": pollutant,
            "value": val,
            "units": units,
            "long_name": long_name,
        })
    return rows


def plot_timeseries_by_pollutant(df_ts: pd.DataFrame, out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)

    for (region, pollutant), g in df_ts.groupby(["region", "pollutant"]):
        g = g.sort_values("year")
        units = g["units"].dropna().iloc[0] if g["units"].notna().any() else ""
        long_name = g["long_name"].dropna().iloc[0] if g["long_name"].notna().any() else pollutant

        fig = plt.figure(figsize=(10.5, 5.5), dpi=200)
        ax = fig.add_subplot(111)

        for (sigla, capital), gg in g.groupby(["sigla", "capital"]):
            ax.plot(gg["year"], gg["value"], marker="o", linewidth=2, label=f"{sigla} ({capital})")

        ax.set_title(f"{region} – Time series capoluoghi – {pollutant}", fontsize=13)
        ax.set_xlabel("Anno")
        ax.set_ylabel(f"{long_name} [{units}]")
        ax.grid(True, which="major", linestyle="--", alpha=0.4)
        ax.legend(ncol=2, fontsize=9, frameon=False)

        out_png = out_dir / f"{region}_{pollutant}_capitals_timeseries.png"
        fig.tight_layout()
        fig.savefig(out_png)
        plt.close(fig)


# =========================
# MONTAGE
# =========================
def plot_montage_years_single_pollutant(
    region_name: str,
    pollutant: str,
    year_to_da_reg: Dict[int, xr.DataArray],
    geom,
    out_png: Path,
    ncols: int = 4,
):
    if not year_to_da_reg:
        raise RuntimeError("Nessun anno disponibile per il montage.")

    years = sorted(year_to_da_reg.keys())
    da0 = year_to_da_reg[years[0]]

    lat_name = da0.attrs["__lat_name__"]
    lon_name = da0.attrs["__lon_name__"]
    units = da0.attrs.get("units", "") or "µg/m3"
    long_name = da0.attrs.get("long_name", da0.attrs.get("standard_name", pollutant)) or pollutant

    minx, miny, maxx, maxy = geom.bounds
    pad_x = (maxx - minx) * 0.03
    pad_y = (maxy - miny) * 0.03
    extent = [minx - pad_x, maxx + pad_x, miny - pad_y, maxy + pad_y]

    all_vals = []
    for y in years:
        arr = np.asarray(year_to_da_reg[y].values)
        arr = arr[np.isfinite(arr)]
        if arr.size:
            all_vals.append(arr)

    if not all_vals:
        raise RuntimeError("Tutti gli anni sono NaN dopo clip.")

    all_vals = np.concatenate(all_vals)
    vmin = float(np.nanpercentile(all_vals, 2))
    vmax = float(np.nanpercentile(all_vals, 98))
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
        vmin = float(np.nanmin(all_vals))
        vmax = float(np.nanmax(all_vals))

    n = len(years)
    nrows = int(np.ceil(n / ncols))

    fig = plt.figure(figsize=(4.2 * ncols, 3.7 * nrows), dpi=DPI)

    axes = []
    mesh = None

    for i, y in enumerate(years):
        ax = fig.add_subplot(nrows, ncols, i + 1, projection=ccrs.PlateCarree())
        axes.append(ax)

        ax.set_extent(extent, crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.COASTLINE.with_scale("10m"), linewidth=0.5)
        ax.add_feature(cfeature.BORDERS.with_scale("10m"), linewidth=0.3)

        xmaj = np.arange(np.floor(extent[0]), np.ceil(extent[1]) + 0.001, 1.0)
        ymaj = np.arange(np.floor(extent[2]), np.ceil(extent[3]) + 0.001, 0.5)
        ax.set_xticks(xmaj, crs=ccrs.PlateCarree())
        ax.set_yticks(ymaj, crs=ccrs.PlateCarree())
        ax.xaxis.set_major_formatter(LongitudeFormatter(number_format=".0f", degree_symbol="°"))
        ax.yaxis.set_major_formatter(LatitudeFormatter(number_format=".1f", degree_symbol="°"))
        ax.tick_params(axis="both", which="major", direction="in", top=True, right=True, length=4, labelsize=8)
        ax.grid(True, which="major", linewidth=0.25, linestyle="--", alpha=0.35)

        da = year_to_da_reg[y]
        lats = da[lat_name].values
        lons = da[lon_name].values
        data = np.ma.masked_invalid(da.values)

        mesh = ax.pcolormesh(
            lons, lats, data,
            transform=ccrs.PlateCarree(),
            shading="auto",
            cmap=DEFAULT_CMAP,
            vmin=vmin, vmax=vmax,
            zorder=1,
        )

        ax.add_geometries([geom], crs=ccrs.PlateCarree(),
                          facecolor="none", edgecolor="black", linewidth=1.1, zorder=5)

        plot_province_capitals(ax, region_name)
        ax.set_title(str(y), fontsize=11)

    for j in range(n, nrows * ncols):
        ax = fig.add_subplot(nrows, ncols, j + 1)
        ax.axis("off")

    fig.suptitle(f"{region_name} – {pollutant} – medie annue {years[0]}–{years[-1]}", fontsize=16, y=0.98)

    fig.subplots_adjust(left=0.04, right=0.98, top=0.93, bottom=0.12, wspace=0.15, hspace=0.22)
    cax = fig.add_axes([0.18, 0.07, 0.64, 0.03])
    cb = fig.colorbar(mesh, cax=cax, orientation="horizontal")
    cb.set_label(f"{long_name} [{units}]", fontsize=11)
    cb.ax.tick_params(labelsize=9, direction="in")

    fig.text(0.01, 0.02, FOOTER_LEFT, ha="left", va="bottom", fontsize=9)
    fig.text(0.99, 0.02, FOOTER_RIGHT, ha="right", va="bottom", fontsize=9)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=DPI)
    plt.close(fig)


# =========================
# GIF
# =========================
def build_gif(frames: List[Path], out_gif: Path, fps: float = GIF_FPS, loop: int = GIF_LOOP) -> None:
    out_gif.parent.mkdir(parents=True, exist_ok=True)
    imgs = []
    for p in frames:
        if p.exists():
            imgs.append(imageio.imread(p))
    if not imgs:
        raise RuntimeError(f"Nessun frame valido per GIF: {out_gif}")
    imageio.mimsave(out_gif, imgs, fps=fps, loop=loop)


# =========================
# MAIN
# =========================
def main():
    print(f"[INFO] DATA_DIR: {DATA_DIR}")
    print(f"[INFO] BOUNDARY_FILE: {BOUNDARY_FILE}")
    print(f"[INFO] OUT_DIR: {OUT_DIR}")
    print(f"[INFO] REGRID_MODE={REGRID_MODE} | REGRID_METHOD={REGRID_METHOD}")

    df_annual, df_monthly = scan_cams_files(DATA_DIR)
    print(f"[INFO] Annual files: {len(df_annual)} | Monthly files: {len(df_monthly)}")

    years = sorted(
        set(df_annual["year"].unique().tolist() if not df_annual.empty else [])
        | set(df_monthly["year"].unique().tolist() if not df_monthly.empty else [])
    )
    pols = sorted(
        set(df_annual["pollutant"].unique().tolist() if not df_annual.empty else [])
        | set(df_monthly["pollutant"].unique().tolist() if not df_monthly.empty else [])
    )

    print(f"[INFO] Anni: {years}")
    print(f"[INFO] Inquinanti: {pols[:30]}{' ...' if len(pols) > 30 else ''}")

    region_geoms = load_region_geoms()

    # griglia reference per pollutant
    ref_grid_by_pollutant: Dict[str, Tuple[xr.DataArray, xr.DataArray]] = {}

    results = []
    frames_by_rp: Dict[Tuple[str, str], List[Tuple[int, Path]]] = {}

    combos = set()
    if not df_annual.empty:
        combos |= set(map(tuple, df_annual[["year", "pollutant"]].itertuples(index=False, name=None)))
    if not df_monthly.empty:
        combos |= set(map(tuple, df_monthly[["year", "pollutant"]].itertuples(index=False, name=None)))

    ts_rows = []

    # ==========================
    # 1) MAPPE + CSV + FRAMES GIF + TIMESERIES
    # ==========================
    for (year, pollutant) in sorted(combos):
        da_year = None

        # 1) annual_mean se esiste
        if not df_annual.empty:
            hit = df_annual[(df_annual["year"] == year) & (df_annual["pollutant"] == pollutant)]
            if not hit.empty:
                path = hit.iloc[0]["path"]
                try:
                    da_year = compute_annual_from_single_nc(path)
                except Exception as e:
                    print(f"[WARN] Annual {year} {pollutant} fallito ({path}): {e}")
                    da_year = None

        # 2) fallback mensili
        if da_year is None:
            if df_monthly.empty:
                print(f"[SKIP] {year} {pollutant}: nessun mensile e nessun annual.")
                continue

            grp = df_monthly[(df_monthly["year"] == year) & (df_monthly["pollutant"] == pollutant)]
            if grp.empty:
                print(f"[SKIP] {year} {pollutant}: nessun mensile.")
                continue

            paths = pick_monthly_paths_for_year_pollutant(grp)
            try:
                da_year = compute_yearly_mean_from_monthlies(paths)
            except Exception as e:
                print(f"[SKIP] {year} {pollutant}: errore compute dai mensili -> {e}")
                continue

        # normalize lon
        lon_name = da_year.attrs["__lon_name__"]
        da_year = normalize_longitudes(da_year, lon_name)

        # (A) imposta (eventualmente aggiorna) la griglia reference per questo pollutant
        maybe_set_reference_grid(pollutant, da_year, ref_grid_by_pollutant, mode=REGRID_MODE)

        # (B) regrid sulla griglia reference (se non è già quella)
        ref_lats, ref_lons = ref_grid_by_pollutant[pollutant]
        da_year_rg = regrid_to_reference(da_year, ref_lats, ref_lons, method=REGRID_METHOD)

        if DEBUG_GRID:
            lat_name = da_year_rg.attrs["__lat_name__"]
            lon_name = da_year_rg.attrs["__lon_name__"]
            print("[GRID]", year, pollutant,
                  "shape=", tuple(da_year_rg.shape),
                  "dlat=", _grid_step_1d(da_year_rg[lat_name].values),
                  "dlon=", _grid_step_1d(da_year_rg[lon_name].values))

        # --- TIME SERIES capoluoghi: usa SEMPRE la griglia regridded (coerenza nel tempo) ---
        for region_name in REGIONS:
            ts_rows.extend(
                build_capitals_timeseries_rows(
                    da_year=da_year_rg,
                    year=year,
                    pollutant=pollutant,
                    region_name=region_name,
                    method="nearest",
                )
            )

        lat_name = da_year_rg.attrs["__lat_name__"]
        lon_name = da_year_rg.attrs["__lon_name__"]
        units = da_year_rg.attrs.get("units", "")
        var_name = da_year_rg.attrs.get("__var_name__", "")
        n_months = da_year_rg.attrs.get("__n_months__", "")

        for region_name, geom in region_geoms.items():
            try:
                da_reg = clip_to_region(da_year_rg, lat_name, lon_name, geom)

                out_png = plot_map(da_reg, region_name, year, pollutant, geom)

                value = regional_weighted_mean(da_reg)
                results.append({
                    "region": region_name,
                    "year": year,
                    "pollutant": pollutant,
                    "value": value,
                    "units": units,
                    "var_name": var_name,
                    "n_months": n_months,
                })

                frames_by_rp.setdefault((region_name, pollutant), []).append((year, out_png))
                print(f"[OK] {year} {pollutant} {region_name} -> {out_png} | mean={value:g} {units}")

            except Exception as e:
                print(f"[SKIP] {year} {pollutant} {region_name}: {repr(e)}")

    # CSV
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = OUT_DIR / "regional_annual_means.csv"
    pd.DataFrame(results).sort_values(["region", "year", "pollutant"]).to_csv(csv_path, index=False)
    print(f"[DONE] CSV salvato: {csv_path}")

    # GIF
    GIF_DIR.mkdir(parents=True, exist_ok=True)
    for (region_name, pollutant), frames in frames_by_rp.items():
        frames_sorted = [p for (y, p) in sorted(frames, key=lambda t: t[0])]
        out_gif = GIF_DIR / f"{region_name}_{pollutant}_annual_mean.gif"
        try:
            build_gif(frames_sorted, out_gif, fps=GIF_FPS, loop=GIF_LOOP)
            print(f"[DONE] GIF: {out_gif}")
        except Exception as e:
            print(f"[SKIP] GIF {region_name} {pollutant}: {e}")

    # TIME SERIES output
    df_ts = pd.DataFrame(ts_rows)
    ts_dir = OUT_DIR / "timeseries_capoluoghi"
    plot_timeseries_by_pollutant(df_ts, ts_dir)
    df_ts.to_csv(ts_dir / "capitals_timeseries_values.csv", index=False)
    print(f"[DONE] Timeseries salvate in: {ts_dir}")

    # ==========================
    # 2) MONTAGE multi-anno
    # ==========================
    montage_root = OUT_DIR / "montages"
    montage_root.mkdir(parents=True, exist_ok=True)

    all_years = years
    all_pollutants = pols

    for region_name, geom in region_geoms.items():
        for pollutant in all_pollutants:
            year_to_da: Dict[int, xr.DataArray] = {}

            # se non ho reference grid per quel pollutant (mai processato sopra), salto
            if pollutant not in ref_grid_by_pollutant:
                continue
            ref_lats, ref_lons = ref_grid_by_pollutant[pollutant]

            for y in all_years:
                da_y = None

                if not df_annual.empty:
                    hit = df_annual[(df_annual["year"] == y) & (df_annual["pollutant"] == pollutant)]
                    if not hit.empty:
                        try:
                            da_y = compute_annual_from_single_nc(hit.iloc[0]["path"])
                        except Exception as e:
                            print(f"[SKIP] montage annual {region_name} {pollutant} {y}: {e}")
                            da_y = None

                if da_y is None and not df_monthly.empty:
                    grp = df_monthly[(df_monthly["year"] == y) & (df_monthly["pollutant"] == pollutant)]
                    if not grp.empty:
                        try:
                            paths = pick_monthly_paths_for_year_pollutant(grp)
                            da_y = compute_yearly_mean_from_monthlies(paths)
                        except Exception as e:
                            print(f"[SKIP] montage monthly {region_name} {pollutant} {y}: {e}")
                            da_y = None

                if da_y is None:
                    continue

                lon_name = da_y.attrs["__lon_name__"]
                da_y = normalize_longitudes(da_y, lon_name)

                # regrid su reference grid
                da_y = regrid_to_reference(da_y, ref_lats, ref_lons, method=REGRID_METHOD)

                lat_name = da_y.attrs["__lat_name__"]
                lon_name = da_y.attrs["__lon_name__"]

                try:
                    da_reg = clip_to_region(da_y, lat_name, lon_name, geom)
                    arr = np.asarray(da_reg.values)
                    if arr.size == 0 or np.sum(np.isfinite(arr)) == 0:
                        continue
                    year_to_da[y] = da_reg
                except Exception as e:
                    print(f"[SKIP] montage clip {region_name} {pollutant} {y}: {e}")

            if len(year_to_da) < 2:
                print(f"[SKIP] MONTAGE {region_name} {pollutant}: troppo pochi anni ({len(year_to_da)})")
                continue

            y0, y1 = min(year_to_da.keys()), max(year_to_da.keys())
            out_png = montage_root / region_name / f"{region_name}_{pollutant}_{y0}-{y1}_montage.png"

            try:
                plot_montage_years_single_pollutant(
                    region_name=region_name,
                    pollutant=pollutant,
                    year_to_da_reg=year_to_da,
                    geom=geom,
                    out_png=out_png,
                    ncols=4,
                )
                print(f"[DONE] MONTAGE: {out_png}")
            except Exception as e:
                print(f"[SKIP] MONTAGE {region_name} {pollutant}: {e}")


if __name__ == "__main__":
    main()
