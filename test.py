from pathlib import Path
import geopandas as gpd

BOUNDARY_FILE = Path("boundaries/regioni_italia.json")

def load_region_geometry(region_name: str):
    gdf = gpd.read_file(BOUNDARY_FILE)
    mask = gdf["name"].astype(str).str.lower().eq(region_name.lower())
    reg = gdf.loc[mask]
    if reg.empty:
        raise ValueError(f"Regione '{region_name}' non trovata.")
    return reg.unary_union

for r in ["Abruzzo", "Puglia"]:
    geom = load_region_geometry(r)
    print(r, geom.geom_type, geom.bounds)
