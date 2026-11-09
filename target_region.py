"""Draw a configurable target region and export its satellite-analysis grid.

Run: python target_region.py
Dependencies: numpy, matplotlib, geopandas (see requirements.txt).
Boundaries: Natural Earth Admin 0, 1:50m. Internet needed on first run only.
The default rectangle is a provisional study area, not an official definition
of Central Europe. It includes all surfaces inside the box, including water.
"""
from pathlib import Path
from urllib.request import urlopen
import hashlib
import json

import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
import matplotlib.patheffects as pe
import numpy as np

# One configuration file controls both grid creation and analysis.
BASE = Path(__file__).resolve().parent
CONFIG = json.loads((BASE / "config.json").read_text())
LON_MIN, LAT_MIN, LON_MAX, LAT_MAX = CONFIG["region_bounds_lon_lat"]
GRID_STEP_DEG = CONFIG["grid_step_deg"]

DATA = BASE / "data"
FIGURES = BASE / "figures"
URL = ("https://naturalearth.s3.amazonaws.com/50m_cultural/"
       "ne_50m_admin_0_countries.zip")


def main():
    if not (-180 <= LON_MIN < LON_MAX <= 180
            and -90 < LAT_MIN < LAT_MAX < 90 and GRID_STEP_DEG > 0):
        raise ValueError("Check the region bounds and positive grid spacing.")
    nx_float = (LON_MAX - LON_MIN) / GRID_STEP_DEG
    ny_float = (LAT_MAX - LAT_MIN) / GRID_STEP_DEG
    if not (np.isclose(nx_float, round(nx_float))
            and np.isclose(ny_float, round(ny_float))):
        raise ValueError("Each region span must be divisible by GRID_STEP_DEG.")
    nx, ny = round(nx_float), round(ny_float)
    if nx < 1 or ny < 1:
        raise ValueError("The region must contain at least one grid cell.")
    DATA.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    archive = DATA / "ne_50m_admin_0_countries.zip"
    if not archive.exists():
        print("Downloading Natural Earth country boundaries...")
        try:
            with urlopen(URL, timeout=60) as response:
                content = response.read()
            temporary = archive.with_suffix(".download")
            temporary.write_bytes(content)
            temporary.replace(archive)
        except OSError as exc:
            raise RuntimeError(
                f"Could not download the map. Try again with internet access, "
                f"or save the Natural Earth ZIP manually to {archive}."
            ) from exc
    countries = gpd.read_file(archive).to_crs("EPSG:4326")

    # Use cell centres, not cell edges. Default settings give 40 x 20 = 800 targets.
    lons = LON_MIN + (np.arange(nx) + 0.5) * GRID_STEP_DEG
    lats = LAT_MIN + (np.arange(ny) + 0.5) * GRID_STEP_DEG
    lon_grid, lat_grid = np.meshgrid(lons, lats)
    lon, lat = lon_grid.ravel(), lat_grid.ravel()

    # Relative cell areas on a sphere: d_lon * (sin(lat_north)-sin(lat_south)).
    # Normalise to sum to one for area-weighted regional summaries later.
    half_step = GRID_STEP_DEG / 2
    weight = np.deg2rad(GRID_STEP_DEG) * (
        np.sin(np.deg2rad(lat + half_step))
        - np.sin(np.deg2rad(lat - half_step))
    )
    weight /= weight.sum()
    targets = np.column_stack((np.arange(1, len(lon) + 1), lon, lat, weight))
    csv_path = DATA / "target_grid.csv"
    np.savetxt(csv_path, targets, delimiter=",", comments="",
               header="target_id,longitude_deg,latitude_deg,area_weight",
               fmt=["%d", "%.6f", "%.6f", "%.12f"])

    from paper_figures import setup, target_map
    config, nearby = setup()
    target_map(config, nearby)
    metadata = {
        "bounds_lon_lat": [LON_MIN, LAT_MIN, LON_MAX, LAT_MAX],
        "grid_spacing_deg": GRID_STEP_DEG, "grid_convention": "cell centres",
        "target_count": len(lon), "crs": "EPSG:4326", "land_only": False,
        "weights": "normalised spherical cell areas, not ellipsoidal areas",
        "boundary_source_url": URL,
        "boundary_zip_sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
    }
    (DATA / "study_settings.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Saved {len(lon)} targets: {csv_path}")
    print(f"Overleaf figure: {FIGURES / 'target_region.pdf'}")
    print(f"Preview: {FIGURES / 'target_region.png'}")


if __name__ == "__main__":
    main()
