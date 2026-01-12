#!/usr/bin/env python3
"""
Adds export of each subcatchment as:
  A) separate GeoTIFF per catchment (binary mask or ID values), OR
  B) one multi-band GeoTIFF (each band = one catchment mask)

Note: multi-band GeoTIFF can get huge if you have many catchments.
Separate files are usually safer.
"""

import os
from pathlib import Path
import numpy as np

os.environ["GDAL_DATA"] = r"C:\Users\MRE\anaconda4\envs\Subcatchment\Library\share\gdal"
os.environ["PROJ_LIB"]  = r"C:\Users\MRE\anaconda4\envs\Subcatchment\Library\share\proj"

import rasterio
import geopandas as gpd
from shapely.geometry import Point
import subprocess


# --------------------------------------------------------------------
# INPUTS
# --------------------------------------------------------------------
LDD_PATH = r"E:\pythonProject\spu\SPU_LDD_0.tif"
UAA_KM2_PATH = r"E:\pythonProject\spu\upstream_area_km2.tif"
THRESH_KM2 = 10

OUTLETS_SHP = r"E:\pythonProject\spu\outlets.shp"
SUBCATCH_RASTER = r"E:\pythonProject\spu\subcatch.tif"
SUBCATCH_SHP = r"E:\pythonProject\spu\subcatch.shp"

WBT_EXE = r"C:\Users\MRE\anaconda4\envs\Subcatchment\Library\bin\whitebox_tools.exe"
ESRI_PNTR_PATH = r"E:\pythonProject\spu\d8_esri_pointer.tif"

# NEW: export options
EXPORT_DIR = r"E:\pythonProject\spu\subcatch_tiffs"  # for separate files
EXPORT_SEPARATE_TIFFS = True
EXPORT_MULTIBAND_TIFF = False  # set True if you want a multi-band file
MULTIBAND_TIFF_PATH = r"E:\pythonProject\spu\subcatch_multiband.tif"
WRITE_MASKS_AS = "mask"  # "mask" (0/1) or "id" (id value, nodata elsewhere)
# --------------------------------------------------------------------

assert os.path.exists(WBT_EXE), f"whitebox_tools.exe not found: {WBT_EXE}"

for p in [OUTLETS_SHP, SUBCATCH_RASTER, SUBCATCH_SHP, ESRI_PNTR_PATH, MULTIBAND_TIFF_PATH]:
    Path(os.path.dirname(p)).mkdir(parents=True, exist_ok=True)
Path(EXPORT_DIR).mkdir(parents=True, exist_ok=True)

# -----------------------------------------------------------
# 1. Load UAA (km²) and LDD
# -----------------------------------------------------------
with rasterio.open(UAA_KM2_PATH) as ds:
    uaa_km2 = ds.read(1)
    profile = ds.profile.copy()
    transform = ds.transform
    uaa_nodata = ds.nodata
    height, width = uaa_km2.shape
    uaa_crs = ds.crs

with rasterio.open(LDD_PATH) as ds_ldd:
    ldd = ds_ldd.read(1)
    ldd_nodata = ds_ldd.nodata

if ldd.shape != uaa_km2.shape:
    raise RuntimeError(f"LDD shape {ldd.shape} != UAA shape {uaa_km2.shape}")
if rasterio.open(LDD_PATH).transform != transform:
    raise RuntimeError("LDD and UAA transforms differ (grids not aligned)")
if rasterio.open(LDD_PATH).crs != uaa_crs:
    raise RuntimeError("LDD and UAA CRS differ")

uaa_mask = np.zeros(uaa_km2.shape, dtype=bool)
if uaa_nodata is not None:
    uaa_mask |= (uaa_km2 == uaa_nodata)
uaa_mask |= ~np.isfinite(uaa_km2)

river = (uaa_km2 >= THRESH_KM2) & (~uaa_mask)

# -----------------------------------------------------------
# 3–5. Find outlets (same as your script)
# -----------------------------------------------------------
ldd_dirs = {
    1: ( 1, -1), 2: ( 1,  0), 3: ( 1,  1),
    4: ( 0, -1), 5: None,      6: ( 0,  1),
    7: (-1, -1), 8: (-1,  0),  9: (-1,  1),
}
rev_code = {(dy, dx): code for code, v in ldd_dirs.items() if v is not None for (dy, dx) in [v]}

ldd_valid = np.isin(ldd, list(ldd_dirs.keys()))
if ldd_nodata is not None:
    ldd_valid &= (ldd != ldd_nodata)
ldd_valid &= (ldd != 255)

outlet_points = []
for row in range(height):
    for col in range(width):
        if not river[row, col]:
            continue

        upstream_count = 0
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                if dy == 0 and dx == 0:
                    continue
                r2 = row + dy
                c2 = col + dx
                if not (0 <= r2 < height and 0 <= c2 < width):
                    continue
                if not river[r2, c2]:
                    continue
                if not ldd_valid[r2, c2]:
                    continue

                needed_code = rev_code.get((-dy, -dx))
                if needed_code is not None and ldd[r2, c2] == needed_code:
                    upstream_count += 1

        if upstream_count == 0 or upstream_count >= 2:
            x, y = rasterio.transform.xy(transform, row, col)
            outlet_points.append(Point(x, y))

gdf = gpd.GeoDataFrame(
    {"id": np.arange(1, len(outlet_points) + 1, dtype=int)},
    geometry=outlet_points,
    crs=profile["crs"],
)
gdf.to_file(OUTLETS_SHP)
print(f"Saved {len(outlet_points)} outlet points -> {OUTLETS_SHP}")

# -----------------------------------------------------------
# 6. Convert LISFLOOD LDD -> ESRI pointer
# -----------------------------------------------------------
LISFLOOD_TO_ESRI = {1: 8, 2: 4, 3: 2, 4: 16, 5: 0, 6: 1, 7: 32, 8: 64, 9: 128}
PNTR_NODATA = np.uint16(65535)

esri = np.zeros(ldd.shape, dtype=np.uint16)
mask_ldd_nodata = np.zeros(ldd.shape, dtype=bool)
if ldd_nodata is not None:
    mask_ldd_nodata |= (ldd == ldd_nodata)
mask_ldd_nodata |= (ldd == 255)

for k, v in LISFLOOD_TO_ESRI.items():
    esri[ldd == k] = np.uint16(v)
esri[mask_ldd_nodata] = PNTR_NODATA

pntr_profile = profile.copy()
pntr_profile.update(dtype="uint16", nodata=int(PNTR_NODATA), count=1)

with rasterio.open(ESRI_PNTR_PATH, "w", **pntr_profile) as dst:
    dst.write(esri, 1)

# -----------------------------------------------------------
# 7. Watershed delineation
# -----------------------------------------------------------
subprocess.run(
    [
        WBT_EXE,
        "--run=Watershed",
        f"--d8_pntr={ESRI_PNTR_PATH}",
        f"--pour_pts={OUTLETS_SHP}",
        f"--output={SUBCATCH_RASTER}",
        "--esri_pntr",
        "-v",
    ],
    check=True
)
print(f"Subcatchment raster written to: {SUBCATCH_RASTER}")

# -----------------------------------------------------------
# 8. Polygon export (optional)
# -----------------------------------------------------------
subprocess.run(
    [
        WBT_EXE,
        "--run=RasterToVectorPolygons",
        f"--input={SUBCATCH_RASTER}",
        f"--output={SUBCATCH_SHP}",
        "-v",
    ],
    check=True
)
print(f"Subcatchment polygons written to: {SUBCATCH_SHP}")

# -----------------------------------------------------------
# 9. NEW: Export each subcatchment to separate TIFFs or one multiband TIFF
# -----------------------------------------------------------
with rasterio.open(SUBCATCH_RASTER) as ds:
    sub = ds.read(1)
    sub_profile = ds.profile.copy()
    sub_nodata = ds.nodata

# Collect unique IDs (exclude nodata/0)
ids = np.unique(sub[np.isfinite(sub)])
if sub_nodata is not None:
    ids = ids[ids != sub_nodata]
ids = ids[ids != 0]
ids = ids.astype(int)

print(f"Found {len(ids)} subcatchments (unique IDs).")

# Output profile for masks
out_profile = sub_profile.copy()
out_profile.update(compress="lzw")

if WRITE_MASKS_AS.lower() == "mask":
    out_dtype = "uint8"
    out_nodata = 0
elif WRITE_MASKS_AS.lower() == "id":
    out_dtype = sub_profile["dtype"]
    out_nodata = sub_nodata if sub_nodata is not None else 0
else:
    raise ValueError("WRITE_MASKS_AS must be 'mask' or 'id'")

# --- A) Separate files ---
if EXPORT_SEPARATE_TIFFS:
    sep_profile = out_profile.copy()
    sep_profile.update(dtype=out_dtype, nodata=out_nodata, count=1)

    for cid in ids:
        out_path = os.path.join(EXPORT_DIR, f"subcatch_{cid:05d}.tif")

        if WRITE_MASKS_AS.lower() == "mask":
            arr = (sub == cid).astype(np.uint8)
        else:
            arr = np.where(sub == cid, sub.astype(sep_profile["dtype"]), out_nodata)

        with rasterio.open(out_path, "w", **sep_profile) as dst:
            dst.write(arr, 1)

    print(f"Wrote {len(ids)} separate GeoTIFFs to: {EXPORT_DIR}")

# --- B) One multiband file ---
if EXPORT_MULTIBAND_TIFF:
    # Warning: can be very large if many subcatchments
    mb_profile = out_profile.copy()
    mb_profile.update(dtype="uint8", nodata=0, count=len(ids))

    with rasterio.open(MULTIBAND_TIFF_PATH, "w", **mb_profile) as dst:
        for band_i, cid in enumerate(ids, start=1):
            dst.write((sub == cid).astype(np.uint8), band_i)

    print(f"Wrote multiband GeoTIFF: {MULTIBAND_TIFF_PATH} (bands={len(ids)})")
