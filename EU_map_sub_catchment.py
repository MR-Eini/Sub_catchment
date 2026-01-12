#!/usr/bin/env python3
"""
EU subcatchments (large raster-safe):
- Inputs (EU domain):
    1) LDD Europe-masked (LISFLOOD keypad codes 1..9, pit=5, nodata often 255)
    2) Upstream area in km² (computed on FULL domain, then masked/cropped to Europe)

- Steps:
    A) stream mask from UAA threshold
    B) pour points = headwaters + confluences (vectorized; no Python pixel loops)
    C) WhiteboxTools Watershed -> subcatchment ID raster
    D) Optional: polygons
    E) Export each subcatchment as separate cropped GeoTIFFs (recommended) OR multiband (not recommended for EU)

Requirements:
  pip install rasterio geopandas shapely
  WhiteboxTools CLI (whitebox_tools.exe)
"""

import os
from pathlib import Path
import numpy as np

# Set BEFORE importing rasterio (Windows/conda)
os.environ["GDAL_DATA"] = r"C:\Users\MRE\anaconda4\envs\Subcatchment\Library\share\gdal"
os.environ["PROJ_LIB"]  = r"C:\Users\MRE\anaconda4\envs\Subcatchment\Library\share\proj"

import rasterio
import geopandas as gpd
import subprocess
from rasterio.features import rasterize


# -------------------------
# INPUTS (EU)
# -------------------------
OUT_DIR = r"E:\pythonProject\spu\EUROPE"

LDD_PATH     = os.path.join(OUT_DIR, "ldd_europe_masked.tif")
UAA_KM2_PATH  = os.path.join(OUT_DIR, "upstream_area_km2_europe.tif")  # must be km² (full-domain acc, EU-masked)

ESRI_PNTR_PATH = os.path.join(OUT_DIR, "d8_esri_pointer_europe.tif")   # already created by your previous script
WBT_EXE = r"C:\Users\MRE\anaconda4\envs\Subcatchment\Library\bin\whitebox_tools.exe"

THRESH_KM2 = 1000.0  # stream definition

# Pour points output (GeoPackage recommended for many points)
OUTLETS_GPKG = os.path.join(OUT_DIR, "outlets.gpkg")
OUTLETS_LAYER = "outlets"

# Watershed outputs
SUBCATCH_RASTER = os.path.join(OUT_DIR, "subcatch.tif")
SUBCATCH_SHP    = os.path.join(OUT_DIR, "subcatch.gpkg")   # use GeoPackage for huge outputs
SUBCATCH_LAYER  = "subcatch"

# Export options
EXPORT_DIR = os.path.join(OUT_DIR, "subcatch_tiffs")
EXPORT_SEPARATE_TIFFS = True
EXPORT_MULTIBAND_TIFF = False  # EU-scale: usually NOT practical
MULTIBAND_TIFF_PATH = os.path.join(OUT_DIR, "subcatch_multiband.tif")

WRITE_MASKS_AS = "mask"  # "mask" => uint8 0/1 ; "id" => keep id values, nodata elsewhere
CROPPED_EXPORT = True    # IMPORTANT for EU: export each catchment cropped to its bbox
MAX_EXPORT = None        # set e.g. 5000 to limit outputs for testing
# -------------------------

assert os.path.exists(WBT_EXE), f"whitebox_tools.exe not found: {WBT_EXE}"
for p in [OUT_DIR, EXPORT_DIR]:
    Path(p).mkdir(parents=True, exist_ok=True)


# -------------------------
# Helpers
# -------------------------
def run_wbt(args, wd):
    wd = os.path.abspath(wd)
    env = dict(os.environ)
    env["RUST_BACKTRACE"] = "1"
    cmd = [WBT_EXE, f"--wd={wd}"] + args + ["-v"]
    p = subprocess.run(cmd, cwd=wd, env=env, text=True, capture_output=True)
    if p.stdout:
        print(p.stdout)
    if p.stderr:
        print(p.stderr)
    if p.returncode != 0:
        raise RuntimeError(f"WhiteboxTools failed (exit={p.returncode}). Command:\n{cmd}")


def compute_inflow_count(ldd, stream, ldd_valid):
    """
    Vectorized inflow count for each cell (how many neighbouring stream cells flow into it).
    LISFLOOD keypad codes:
      1 SW,2 S,3 SE,4 W,5 pit,6 E,7 NW,8 N,9 NE
    A neighbor flows into current if neighbor's code points from neighbor -> current.
    """
    h, w = ldd.shape

    # downstream offset -> code
    ldd_dirs = {
        1: ( 1, -1), 2: ( 1,  0), 3: ( 1,  1),
        4: ( 0, -1), 5: None,      6: ( 0,  1),
        7: (-1, -1), 8: (-1,  0),  9: (-1,  1),
    }
    rev_code = {(dy, dx): code for code, off in ldd_dirs.items() if off is not None for (dy, dx) in [off]}

    inflow = np.zeros((h, w), dtype=np.uint8)

    for ndy in (-1, 0, 1):
        for ndx in (-1, 0, 1):
            if ndy == 0 and ndx == 0:
                continue

            needed = rev_code.get((-ndy, -ndx))  # neighbor must point to current
            if needed is None:
                continue

            # current indices where neighbor exists
            r0 = max(0, -ndy)
            r1 = h - max(0, ndy)
            c0 = max(0, -ndx)
            c1 = w - max(0, ndx)

            dst = (slice(r0, r1), slice(c0, c1))
            src = (slice(r0 + ndy, r1 + ndy), slice(c0 + ndx, c1 + ndx))

            nbr_ok = (
                stream[src]
                & ldd_valid[src]
                & (ldd[src] == needed)
            )

            inflow[dst] += nbr_ok.astype(np.uint8)

    return inflow


def save_points(rows, cols, transform, crs, out_gpkg, layer):
    # affine transform (handles north-up grids)
    xs = transform.c + (cols + 0.5) * transform.a + (rows + 0.5) * transform.b
    ys = transform.f + (cols + 0.5) * transform.d + (rows + 0.5) * transform.e

    gdf = gpd.GeoDataFrame(
        {"id": np.arange(1, len(rows) + 1, dtype=np.int64)},
        geometry=gpd.points_from_xy(xs, ys),
        crs=crs,
    )
    # overwrite layer if exists
    if os.path.exists(out_gpkg):
        # geopandas will append layers; to overwrite reliably, delete file for simplicity
        os.remove(out_gpkg)
    gdf.to_file(out_gpkg, layer=layer, driver="GPKG")
    return len(gdf)

def export_catchments_cropped(subcatch_tif, export_dir, write_as="mask", max_export=None):
    Path(export_dir).mkdir(parents=True, exist_ok=True)

    with rasterio.open(subcatch_tif) as ds:
        sub = ds.read(1)
        prof0 = ds.profile.copy()
        tfm0 = ds.transform
        nd = ds.nodata

    valid = np.isfinite(sub)
    if nd is not None:
        valid &= (sub != nd)
    valid &= (sub != 0)

    rows, cols = np.nonzero(valid)
    vals = sub[rows, cols].astype(np.int64)

    ids, inv = np.unique(vals, return_inverse=True)

    if max_export is not None and len(ids) > max_export:
        ids = ids[:max_export]
        keep = np.isin(vals, ids)
        rows, cols, vals = rows[keep], cols[keep], vals[keep]
        ids, inv = np.unique(vals, return_inverse=True)

    print(f"Exporting {len(ids)} catchments to: {export_dir}")

    n = len(ids)
    min_r = np.full(n, np.iinfo(np.int32).max, dtype=np.int32)
    max_r = np.full(n, -1, dtype=np.int32)
    min_c = np.full(n, np.iinfo(np.int32).max, dtype=np.int32)
    max_c = np.full(n, -1, dtype=np.int32)

    np.minimum.at(min_r, inv, rows.astype(np.int32))
    np.maximum.at(max_r, inv, rows.astype(np.int32))
    np.minimum.at(min_c, inv, cols.astype(np.int32))
    np.maximum.at(max_c, inv, cols.astype(np.int32))

    if write_as.lower() == "mask":
        out_dtype = "uint8"
        out_nodata = 0
    elif write_as.lower() == "id":
        out_dtype = prof0["dtype"]
        out_nodata = nd if nd is not None else 0
    else:
        raise ValueError("write_as must be 'mask' or 'id'")

    def tile_dim(dim, preferred=512):
        """Return a valid tile size (multiple of 16) <= dim, or None if tiling not possible."""
        if dim < 16:
            return None
        t = min(preferred, dim)
        t = (t // 16) * 16
        return max(16, t)

    for i, cid in enumerate(ids):
        r0, r1 = int(min_r[i]), int(max_r[i]) + 1
        c0, c1 = int(min_c[i]), int(max_c[i]) + 1
        if r0 >= r1 or c0 >= c1:
            continue

        window = rasterio.windows.Window.from_slices((r0, r1), (c0, c1))
        out_tfm = rasterio.windows.transform(window, tfm0)

        if write_as.lower() == "mask":
            arr = (sub[r0:r1, c0:c1] == cid).astype(np.uint8)
        else:
            arr = np.where(sub[r0:r1, c0:c1] == cid, sub[r0:r1, c0:c1], out_nodata)

        h, w = arr.shape
        bx = tile_dim(w)
        by = tile_dim(h)

        out_prof = prof0.copy()
        out_prof.update(
            driver="GTiff",
            height=h,
            width=w,
            transform=out_tfm,
            dtype=out_dtype,
            nodata=out_nodata,
            count=1,
            compress="lzw",
            BIGTIFF="IF_SAFER",
        )

        # Only tile when we can use valid multiples-of-16 tiles
        if bx is not None and by is not None:
            out_prof.update(tiled=True, blockxsize=bx, blockysize=by)
        else:
            out_prof.update(tiled=False)
            # Ensure no stale block sizes remain
            out_prof.pop("blockxsize", None)
            out_prof.pop("blockysize", None)

        out_path = os.path.join(export_dir, f"subcatch_{int(cid):08d}.tif")
        with rasterio.open(out_path, "w", **out_prof) as dst:
            dst.write(arr.astype(out_dtype), 1)

        if (i + 1) % 500 == 0:
            print(f"  ... {i+1}/{len(ids)} written")

    print("Done exporting catchments.")


# -------------------------
# Main
# -------------------------
print("Reading EU LDD + UAA...")
with rasterio.open(UAA_KM2_PATH) as ds:
    uaa_km2 = ds.read(1)
    transform = ds.transform
    crs = ds.crs
    uaa_nd = ds.nodata
    prof_u = ds.profile.copy()

with rasterio.open(LDD_PATH) as ds:
    ldd = ds.read(1)
    ldd_nd = ds.nodata
    prof_l = ds.profile.copy()

# Alignment checks
if ldd.shape != uaa_km2.shape:
    raise RuntimeError(f"Shape mismatch: LDD {ldd.shape} vs UAA {uaa_km2.shape}")
if prof_l["transform"] != prof_u["transform"]:
    raise RuntimeError("Transform mismatch: LDD and UAA grids not aligned")
if prof_l.get("crs") != prof_u.get("crs"):
    raise RuntimeError("CRS mismatch: LDD and UAA differ")

# Masks
uaa_mask = ~np.isfinite(uaa_km2)
if uaa_nd is not None:
    uaa_mask |= (uaa_km2 == uaa_nd)

ldd_valid = np.isin(ldd, [1,2,3,4,5,6,7,8,9])
if ldd_nd is not None:
    ldd_valid &= (ldd != ldd_nd)
ldd_valid &= (ldd != 255)

# Stream raster from upstream area threshold
stream = (uaa_km2 >= float(THRESH_KM2)) & (~uaa_mask) & ldd_valid

print("Computing inflow counts (vectorized)...")
inflow = compute_inflow_count(ldd, stream, ldd_valid)

# Headwaters (0 inflow) or confluences (>=2 inflow) on stream cells
pour_mask = stream & ((inflow == 0) | (inflow >= 2))
rows, cols = np.nonzero(pour_mask)

print(f"Pour points: {len(rows)}")
print("Writing pour points (GeoPackage)...")
npts = save_points(rows, cols, transform, crs, OUTLETS_GPKG, OUTLETS_LAYER)
print(f"Wrote {npts} points -> {OUTLETS_GPKG} (layer={OUTLETS_LAYER})")

POUR_RASTER = os.path.join(OUT_DIR, "pour_pts.tif")

# Read the outlets you just wrote (or keep gdf from save_points if you return it)
gdf = gpd.read_file(OUTLETS_GPKG, layer=OUTLETS_LAYER)

with rasterio.open(ESRI_PNTR_PATH) as tpl:
    out_shape = (tpl.height, tpl.width)
    tfm = tpl.transform
    crs = tpl.crs
    prof = tpl.profile.copy()

# Rasterize: value = outlet id, background = 0
shapes = ((geom, int(fid)) for geom, fid in zip(gdf.geometry, gdf["id"]))
pour = rasterize(
    shapes=shapes,
    out_shape=out_shape,
    transform=tfm,
    fill=0,
    dtype="int32"
)

prof.update(
    driver="GTiff",
    dtype="int32",
    count=1,
    nodata=0,
    crs=crs,
    transform=tfm,
    compress="lzw",
    tiled=True,
    blockxsize=min(512, out_shape[1]),
    blockysize=min(512, out_shape[0]),
    BIGTIFF="IF_SAFER",
)

with rasterio.open(POUR_RASTER, "w", **prof) as dst:
    dst.write(pour, 1)

print("Wrote pour-point raster:", POUR_RASTER)


# Watershed delineation
print("Running Whitebox Watershed...")
run_wbt(
    [
        "--run=Watershed",
        f"--d8_pntr={ESRI_PNTR_PATH}",
        f"--pour_pts={POUR_RASTER}",      # Whitebox reads raster fine
        f"--output={SUBCATCH_RASTER}",
        "--esri_pntr",
    ],
    wd=OUT_DIR
)
print("Wrote subcatchment raster:", SUBCATCH_RASTER)

# Optional polygons (can be very large)
print("Raster -> polygons (optional, may be heavy on EU)...")
run_wbt(
    [
        "--run=RasterToVectorPolygons",
        f"--input={SUBCATCH_RASTER}",
        f"--output={SUBCATCH_SHP}",
    ],
    wd=OUT_DIR
)
print("Wrote subcatchment polygons:", SUBCATCH_SHP)

# Export each catchment
if EXPORT_SEPARATE_TIFFS:
    if not CROPPED_EXPORT:
        raise RuntimeError("For EU-scale you should keep CROPPED_EXPORT=True (full-grid per ID is not feasible).")
    export_catchments_cropped(
        SUBCATCH_RASTER,
        EXPORT_DIR,
        write_as=WRITE_MASKS_AS,
        max_export=MAX_EXPORT
    )

# Multiband export (generally impractical for EU-scale, but included)
if EXPORT_MULTIBAND_TIFF:
    with rasterio.open(SUBCATCH_RASTER) as ds:
        sub = ds.read(1)
        sp = ds.profile.copy()
        nd = ds.nodata

    ids = np.unique(sub[np.isfinite(sub)])
    if nd is not None:
        ids = ids[ids != nd]
    ids = ids[ids != 0].astype(int)

    if MAX_EXPORT is not None and len(ids) > MAX_EXPORT:
        ids = ids[:MAX_EXPORT]

    mb_profile = sp.copy()
    mb_profile.update(driver="GTiff", dtype="uint8", nodata=0, count=len(ids), compress="lzw", BIGTIFF="IF_SAFER")

    with rasterio.open(MULTIBAND_TIFF_PATH, "w", **mb_profile) as dst:
        for band_i, cid in enumerate(ids, start=1):
            dst.write((sub == cid).astype(np.uint8), band_i)

    print(f"Wrote multiband GeoTIFF: {MULTIBAND_TIFF_PATH} (bands={len(ids)})")
