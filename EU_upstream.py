#!/usr/bin/env python3
"""
Compute upstream area (km²) from a very large LISFLOOD/PCRaster LDD .map by:
  1) masking/cropping to Europe (polygon-based; bbox fallback),
  2) converting LISFLOOD keypad LDD -> ESRI D8 pointer (powers of two),
  3) WhiteboxTools D8FlowAccumulation (out_type=cells),
  4) converting upstream cells -> km² (EPSG:4326, latitude-dependent cell area).

Outputs:
  - ldd_europe_masked.tif
  - d8_esri_pointer_europe.tif
  - upstream_cells_europe.tif
  - upstream_area_km2_europe.tif
"""

import os
import math
from pathlib import Path
import numpy as np

# Set BEFORE importing rasterio (helps with GDAL/PROJ issues on Windows/conda)
os.environ["GDAL_DATA"] = r"C:\Users\MRE\anaconda4\envs\Subcatchment\Library\share\gdal"
os.environ["PROJ_LIB"]  = r"C:\Users\MRE\anaconda4\envs\Subcatchment\Library\share\proj"

import rasterio
from rasterio.features import geometry_mask
from rasterio.windows import from_bounds
import subprocess

# Optional but recommended for a true "Europe" mask polygon
try:
    import geopandas as gpd
    from shapely.geometry import box
except Exception:
    gpd = None
    box = None


# --------------------------------------------------------------------
# INPUTS
# --------------------------------------------------------------------
LDD_MAP_PATH = r"E:\pythonProject\spu\LDD_EU_01Min.tif"
WBT_EXE = r"C:\Users\MRE\anaconda4\envs\Subcatchment\Library\bin\whitebox_tools.exe"
assert os.path.exists(WBT_EXE), f"whitebox_tools.exe not found: {WBT_EXE}"

OUT_DIR = r"E:\pythonProject\spu\EUROPE"
Path(OUT_DIR).mkdir(parents=True, exist_ok=True)

LDD_EU_TIF       = os.path.join(OUT_DIR, "ldd_europe_masked.tif")
ESRI_PNTR_EU_TIF = os.path.join(OUT_DIR, "d8_esri_pointer_europe.tif")
UP_CELLS_EU_TIF  = os.path.join(OUT_DIR, "upstream_cells_europe.tif")
UAA_KM2_EU_TIF   = os.path.join(OUT_DIR, "upstream_area_km2_europe.tif")

# Europe bbox (edit if needed)
EU_BBOX = (-25.0, 34.0, 45.0, 72.0)  # (min_lon, min_lat, max_lon, max_lat)
# --------------------------------------------------------------------


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

# keypad code -> downstream offset (dy, dx) in row/col
LDD_DIRS = {
    1: ( 1, -1),  # SW
    2: ( 1,  0),  # S
    3: ( 1,  1),  # SE
    4: ( 0, -1),  # W
    5: None,      # pit/outlet
    6: ( 0,  1),  # E
    7: (-1, -1),  # NW
    8: (-1,  0),  # N
    9: (-1,  1),  # NE
}

PNTR_NODATA = np.uint16(65535)


def get_europe_geometry():
    """Return a Europe polygon (shapely geometry). Uses Natural Earth lowres if available; bbox fallback."""
    minx, miny, maxx, maxy = EU_BBOX

    # Always apply bbox in the end (prevents overseas territories, reduces extent)
    if box is None:
        # no shapely available -> bbox as a tuple-like polygon not possible; use bbox-only logic later
        return None

    europe_bbox_geom = box(minx, miny, maxx, maxy)

    if gpd is None:
        return europe_bbox_geom  # fallback

    try:
        # GeoPandas built-in dataset (often available)
        world_path = gpd.datasets.get_path("naturalearth_lowres")
        world = gpd.read_file(world_path)

        # union all countries labeled as Europe
        eu = world[world["continent"] == "Europe"].to_crs("EPSG:4326")
        geom = eu.unary_union

        # clip to bbox to remove overseas territories (e.g., French Guiana)
        geom = geom.intersection(europe_bbox_geom)
        return geom
    except Exception:
        return europe_bbox_geom  # fallback


def write_masked_europe_ldd(src_path: str, out_tif: str) -> None:
    """Crop to Europe bbox and mask to Europe polygon, writing a GeoTIFF LDD."""
    minx, miny, maxx, maxy = EU_BBOX
    europe_geom = get_europe_geometry()

    with rasterio.open(src_path) as src:
        crs = src.crs
        if crs is None:
            # Many PCRaster .map files are geographic; assume EPSG:4326 if missing
            crs = rasterio.crs.CRS.from_epsg(4326)

        nodata = src.nodata
        if nodata is None:
            # LISFLOOD commonly uses 255 as NoData for LDD
            nodata = 255

        # Compute a window for bbox, read only that window (fast, avoids full raster read)
        win = from_bounds(minx, miny, maxx, maxy, transform=src.transform)
        win = win.round_offsets().round_lengths()

        ldd = src.read(1, window=win)
        tfm = src.window_transform(win)

        # Mask to Europe geometry (polygon if available; bbox otherwise)
        if europe_geom is not None:
            inside = geometry_mask([europe_geom], out_shape=ldd.shape, transform=tfm, invert=True)
            ldd = np.where(inside, ldd, nodata).astype(ldd.dtype)

        profile = src.profile.copy()
        profile.update(
            driver="GTiff",
            height=ldd.shape[0],
            width=ldd.shape[1],
            transform=tfm,
            crs=crs,
            nodata=nodata,
            compress="lzw",
            tiled=True,
            blockxsize=512,
            blockysize=512,
            count=1,
        )

        with rasterio.open(out_tif, "w", **profile) as dst:
            dst.write(ldd, 1)


def build_esri_pointer(ldd_tif: str, out_pntr_tif: str) -> None:
    """
    Convert LISFLOOD keypad LDD -> ESRI D8 pointer GeoTIFF and VERIFY it exists.
    This version writes using a minimal explicit GTiff profile and uses a temp file + rename.
    """
    out_path = Path(out_pntr_tif)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with rasterio.open(ldd_tif) as src:
        ldd = src.read(1)
        tfm = src.transform
        crs = src.crs
        ldd_nodata = src.nodata

    if ldd_nodata is None:
        ldd_nodata = 255

    h, w = ldd.shape

    # inside-domain cells: valid LDD codes and not nodata
    valid_codes = np.array(list(LISFLOOD_TO_ESRI.keys()), dtype=ldd.dtype)
    inside = np.isin(ldd, valid_codes) & (ldd != ldd_nodata) & (ldd != 255)

    esri = np.zeros((h, w), dtype=np.uint16)
    for k, v in LISFLOOD_TO_ESRI.items():
        esri[ldd == k] = np.uint16(v)

    # pointer nodata must be distinct from 0
    esri[~inside] = PNTR_NODATA

    # Write temp then rename
    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")
    if tmp_path.exists():
        tmp_path.unlink()
    if out_path.exists():
        out_path.unlink()

    blockx = min(512, w)
    blocky = min(512, h)

    print("   Writing pointer temp:", str(tmp_path))

    with rasterio.open(
        tmp_path,
        "w",
        driver="GTiff",
        height=h,
        width=w,
        count=1,
        dtype="uint16",
        crs=crs,
        transform=tfm,
        nodata=int(PNTR_NODATA),
        compress="lzw",
        tiled=True,
        blockxsize=blockx,
        blockysize=blocky,
        BIGTIFF="IF_SAFER",
    ) as dst:
        dst.write(esri, 1)

    tmp_path.replace(out_path)

    # Hard verification
    if not out_path.exists() or out_path.stat().st_size == 0:
        raise RuntimeError(f"Pointer raster not created: {out_pntr_tif}")

    print("   Pointer written OK. Size (bytes):", out_path.stat().st_size)


def run_wbt(tool_args, wd):
    wd = os.path.abspath(wd)
    assert os.path.isdir(wd), f"Whitebox working dir does not exist: {wd}"

    env = dict(os.environ)
    env["RUST_BACKTRACE"] = "1"   # useful if Whitebox panics again

    cmd = [WBT_EXE, f"--wd={wd}"] + tool_args + ["-v"]  # --wd supported by CLI :contentReference[oaicite:1]{index=1}

    p = subprocess.run(
        cmd,
        cwd=wd,
        env=env,
        text=True,
        capture_output=True
    )

    # Always print output (so you see what happened in PyCharm)
    if p.stdout:
        print(p.stdout)
    if p.stderr:
        print(p.stderr)

    if p.returncode != 0:
        raise RuntimeError(f"WhiteboxTools failed (exit={p.returncode}). Command:\n{cmd}")

def flow_accum_cells(pntr_tif: str, out_cells_tif: str) -> None:
    assert os.path.exists(pntr_tif), f"Pointer raster not found: {pntr_tif}"
    Path(os.path.dirname(out_cells_tif)).mkdir(parents=True, exist_ok=True)

    run_wbt(
        [
            "--run=D8FlowAccumulation",
            f"--input={pntr_tif}",
            f"--output={out_cells_tif}",
            "--out_type=cells",
            "--pntr",
            "--esri_pntr",
        ],
        wd=OUT_DIR
    )

def cells_to_km2(up_cells_tif: str, out_km2_tif: str) -> None:
    """Convert upstream cell counts to upstream area in km² for EPSG:4326 rasters."""
    with rasterio.open(up_cells_tif) as ds:
        up = ds.read(1).astype("float64")
        prof = ds.profile.copy()
        nd = ds.nodata
        tfm = ds.transform
        crs = ds.crs
        h, w = up.shape

        if crs is None or not crs.is_geographic:
            raise RuntimeError(f"Expected geographic CRS (EPSG:4326-like), got: {crs}")

        dlon_deg = float(tfm.a)
        dlat_deg = abs(float(tfm.e))

        rows = np.arange(h, dtype="float64")
        lat_deg = tfm.f + (rows + 0.5) * tfm.e
        lat_rad = np.deg2rad(lat_deg)

        R = 6378137.0
        dlon_rad = math.radians(dlon_deg)
        dlat_rad = math.radians(dlat_deg)

        cell_km2_per_row = (R * R * dlat_rad * dlon_rad * np.cos(lat_rad)) / 1e6
        uaa_km2 = up * cell_km2_per_row[:, None]

        if nd is not None:
            uaa_km2[up == nd] = nd

        prof.update(dtype="float32", compress="lzw", tiled=True, blockxsize=512, blockysize=512)

        with rasterio.open(out_km2_tif, "w", **prof) as out:
            out.write(uaa_km2.astype("float32"), 1)

import time

def mask_raster_to_europe(src_tif: str, out_tif: str) -> None:
    minx, miny, maxx, maxy = EU_BBOX
    europe_geom = get_europe_geometry()

    src_abs = os.path.abspath(src_tif)
    out_abs = os.path.abspath(out_tif)
    out_path = Path(out_abs)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Always write to temp first (prevents overwrite issues)
    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")

    # --- Read first (and CLOSE) ---
    with rasterio.open(src_abs) as src:
        win = from_bounds(minx, miny, maxx, maxy, transform=src.transform)
        win = win.round_offsets().round_lengths()

        arr = src.read(1, window=win)
        tfm = src.window_transform(win)
        profile = src.profile.copy()
        nodata = src.nodata

        if nodata is None:
            nodata = -9999.0

        if europe_geom is not None:
            inside = geometry_mask([europe_geom], out_shape=arr.shape, transform=tfm, invert=True)
            arr = np.where(inside, arr, nodata)

    # Update profile for output window
    profile.update(
        driver="GTiff",
        height=arr.shape[0],
        width=arr.shape[1],
        transform=tfm,
        nodata=nodata,
        compress="lzw",
        tiled=True,
        blockxsize=min(512, arr.shape[1]),
        blockysize=min(512, arr.shape[0]),
        count=1,
        BIGTIFF="IF_SAFER",
    )

    # Remove temp if exists
    if tmp_path.exists():
        tmp_path.unlink()

    # --- Write temp ---
    with rasterio.open(tmp_path, "w", **profile) as dst:
        dst.write(arr.astype(profile["dtype"]), 1)

    # --- Replace final ---
    # If output exists and is locked by another program, you'll still get PermissionError here.
    for _ in range(5):
        try:
            if out_path.exists():
                out_path.unlink()
            tmp_path.replace(out_path)
            break
        except PermissionError:
            time.sleep(0.5)
    else:
        raise PermissionError(
            f"Cannot overwrite output (locked). Close QGIS/ArcGIS/PyCharm viewer and retry:\n{out_abs}"
        )



# ---------------------------
# Run
# ---------------------------
print("1) Convert FULL LDD to ESRI D8 pointer...")
build_esri_pointer(LDD_MAP_PATH, ESRI_PNTR_EU_TIF)  # pointer for full domain OR write to a full-domain path
print("   Wrote:", ESRI_PNTR_EU_TIF)

print("2) Flow accumulation (cells) on FULL domain...")
flow_accum_cells(ESRI_PNTR_EU_TIF, UP_CELLS_EU_TIF)
print("   Wrote:", UP_CELLS_EU_TIF)

print("3) Convert cells -> upstream area (km²) on FULL domain...")
cells_to_km2(UP_CELLS_EU_TIF, UAA_KM2_EU_TIF)
print("   Wrote:", UAA_KM2_EU_TIF)

print("4) Mask/crop Europe for outputs...")
write_masked_europe_ldd(LDD_MAP_PATH, LDD_EU_TIF)  # for LDD Europe output
print("   Wrote:", LDD_EU_TIF)

# Crop/mask upstream_area_km2 to Europe (new helper needed)
mask_raster_to_europe(UAA_KM2_EU_TIF, os.path.join(OUT_DIR, "upstream_area_km2_europe.tif"))
print("   Wrote:", os.path.join(OUT_DIR, "upstream_area_km2_europe.tif"))

print("5) Sanity check...")
with rasterio.open(UAA_KM2_EU_TIF) as ds:
    a = ds.read(1)
    nd = ds.nodata
    print("   CRS:", ds.crs)
    print("   RES:", ds.res)
    if nd is not None:
        a = a[a != nd]
    a = a[np.isfinite(a)]
    print("   min km2:", float(a.min()))
    print("   max km2:", float(a.max()))
