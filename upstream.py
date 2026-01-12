import os
import math
from pathlib import Path
import numpy as np

# Set BEFORE importing rasterio in some setups
os.environ["GDAL_DATA"] = r"C:\Users\MRE\anaconda4\envs\Subcatchment\Library\share\gdal"
os.environ["PROJ_LIB"]  = r"C:\Users\MRE\anaconda4\envs\Subcatchment\Library\share\proj"

import rasterio
import subprocess

LDD_PATH = r"E:\pythonProject\spu\SPU_LDD_0.tif"
ESRI_PNTR_PATH = r"E:\pythonProject\spu\d8_esri_pointer.tif"

UP_CELLS_PATH = r"E:\pythonProject\spu\upstream_cells.tif"
UAA_KM2_PATH  = r"E:\pythonProject\spu\upstream_area_km2.tif"

WBT_EXE = r"C:\Users\MRE\anaconda4\envs\Subcatchment\Library\bin\whitebox_tools.exe"
assert os.path.exists(WBT_EXE), f"whitebox_tools.exe not found: {WBT_EXE}"

# LISFLOOD/PCRaster keypad LDD -> ESRI D8 (powers of two)
LISFLOOD_TO_ESRI = {
    1: 8,    # SW
    2: 4,    # S
    3: 2,    # SE
    4: 16,   # W
    5: 0,    # pit/outlet
    6: 1,    # E
    7: 32,   # NW
    8: 64,   # N
    9: 128,  # NE
}

Path(os.path.dirname(ESRI_PNTR_PATH)).mkdir(parents=True, exist_ok=True)
Path(os.path.dirname(UP_CELLS_PATH)).mkdir(parents=True, exist_ok=True)
Path(os.path.dirname(UAA_KM2_PATH)).mkdir(parents=True, exist_ok=True)

# --- Read LDD and build ESRI pointer raster ---
with rasterio.open(LDD_PATH) as src:
    ldd = src.read(1)
    prof = src.profile.copy()
    ldd_nodata = src.nodata
    transform = src.transform
    crs = src.crs
    h, w = ldd.shape

print("LDD CRS:", crs)
print("LDD RES:", src.res)

esri = np.zeros(ldd.shape, dtype=np.uint16)

PNTR_NODATA = np.uint16(65535)
mask_nodata = np.zeros(ldd.shape, dtype=bool)
if ldd_nodata is not None:
    mask_nodata |= (ldd == ldd_nodata)
mask_nodata |= (ldd == 255)  # common LISFLOOD NoData

for k, v in LISFLOOD_TO_ESRI.items():
    esri[ldd == k] = np.uint16(v)

esri[mask_nodata] = PNTR_NODATA
prof.update(dtype="uint16", nodata=int(PNTR_NODATA), count=1)

with rasterio.open(ESRI_PNTR_PATH, "w", **prof) as dst:
    dst.write(esri, 1)

# --- Flow accumulation as upstream cell counts (robust, no units) ---
subprocess.run(
    [
        WBT_EXE,
        "--run=D8FlowAccumulation",
        f"--input={ESRI_PNTR_PATH}",
        f"--output={UP_CELLS_PATH}",
        "--out_type=cells",
        "--pntr",
        "--esri_pntr",
        "-v",
    ],
    check=True
)

print("Wrote upstream cell-count raster:", UP_CELLS_PATH)

# --- Convert upstream cells -> km² for EPSG:4326 ---
with rasterio.open(UP_CELLS_PATH) as ds:
    up_cells = ds.read(1).astype("float64")
    out_prof = ds.profile.copy()
    nd = ds.nodata
    tfm = ds.transform
    h, w = up_cells.shape

    if ds.crs is None or not ds.crs.is_geographic:
        raise RuntimeError(
            f"Expected geographic CRS (EPSG:4326-like) for km² conversion, got: {ds.crs}"
        )

    # Pixel size in degrees
    dlon_deg = float(tfm.a)
    dlat_deg = abs(float(tfm.e))

    # Latitude at each row center
    rows = np.arange(h, dtype="float64")
    lat_deg = tfm.f + (rows + 0.5) * tfm.e
    lat_rad = np.deg2rad(lat_deg)

    # Cell area in km² per row (WGS84 sphere approximation)
    R = 6378137.0  # meters
    dlon_rad = math.radians(dlon_deg)
    dlat_rad = math.radians(dlat_deg)
    cell_km2_per_row = (R * R * dlat_rad * dlon_rad * np.cos(lat_rad)) / 1e6  # (h,)

    uaa_km2 = up_cells * cell_km2_per_row[:, None]

    if nd is not None:
        uaa_km2[up_cells == nd] = nd

    out_prof.update(dtype="float32")
    with rasterio.open(UAA_KM2_PATH, "w", **out_prof) as out:
        out.write(uaa_km2.astype("float32"), 1)

print("Wrote upstream area (km²) raster:", UAA_KM2_PATH)

# --- Sanity check ---
with rasterio.open(UAA_KM2_PATH) as ds:
    a = ds.read(1)
    nd = ds.nodata
    print("UAA CRS:", ds.crs)
    print("UAA RES:", ds.res)
    if nd is not None:
        a = a[a != nd]
    a = a[np.isfinite(a)]
    print("min km2:", float(a.min()))
    print("max km2:", float(a.max()))
