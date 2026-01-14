#!/usr/bin/env python3
"""
EU subcatchments from DEM (large raster-safe), analogous to your LDD/UAA workflow:

Steps:
  1) Fill depressions (optional but recommended)
  2) D8 pointer from DEM (ESRI-style powers of two)
  3) D8 flow accumulation (cells) on full domain
  4) Convert cells -> upstream area km² (constant cell area for projected CRS)
  5) Stream mask: UAA_km2 >= THRESH_KM2
  6) Pour points raster: headwaters (inflow=0) + confluences (inflow>=2) on stream network (blockwise)
  7) Whitebox Watershed (pointer + pour raster) -> subcatch raster
  8) Optional polygons
  9) Export each catchment as cropped GeoTIFFs (streamed; avoids full reads)

Requirements:
  pip install rasterio numpy
  (optional) geopandas if you want point layers; not required here
  WhiteboxTools CLI (whitebox_tools.exe)
"""

import os
import math
from pathlib import Path
import numpy as np

# If you had GDAL warnings in conda, keep these
os.environ["GDAL_DATA"] = r"C:\Users\MRE\anaconda4\envs\Subcatchment\Library\share\gdal"
os.environ["PROJ_LIB"]  = r"C:\Users\MRE\anaconda4\envs\Subcatchment\Library\share\proj"

import rasterio
from rasterio.windows import Window
import subprocess


# -------------------------
# INPUTS
# -------------------------
OUT_DIR = r"Z:\EUData\download\Mamad\Europe\SubCatch_from_DTM30m\OUT_DEM90M_WHITEBOX"
Path(OUT_DIR).mkdir(parents=True, exist_ok=True)

DEM_90M = r"Z:\EUData\download\Mamad\Europe\SubCatch_from_DTM30m\dem_90m_3035.tif"  # <-- your 90m DEM (create in ArcGIS if needed)

# Optional: use a pre-masked DEM (sea/lakes as NoData) if you have it
DEM_MASKED = None  # e.g. r"...\dem_90m_land_only.tif"

# Optional: use a pre-burned DEM (burn rivers in ArcGIS if needed)
DEM_BURNED = None  # e.g. r"...\dem_90m_burned.tif"

# WhiteboxTools
WBT_EXE = r"C:\Users\MRE\anaconda4\envs\Subcatchment\Library\bin\whitebox_tools.exe"
assert os.path.exists(WBT_EXE), f"whitebox_tools.exe not found: {WBT_EXE}"

THRESH_KM2 = 1000.0  # stream definition threshold (like your LDD/UAA script)

# Outputs / intermediates
DEM_FILLED   = os.path.join(OUT_DIR, "dem_filled.tif")
D8_PNTR      = os.path.join(OUT_DIR, "d8_pointer_esri.tif")
UP_CELLS     = os.path.join(OUT_DIR, "upstream_cells.tif")
UAA_KM2      = os.path.join(OUT_DIR, "upstream_area_km2.tif")
POUR_RASTER  = os.path.join(OUT_DIR, "pour_pts.tif")

SUBCATCH_RASTER = os.path.join(OUT_DIR, "subcatch.tif")
SUBCATCH_GPKG   = os.path.join(OUT_DIR, "subcatch.gpkg")   # optional polygons
SUBCATCH_LAYER  = "subcatch"

EXPORT_DIR = os.path.join(OUT_DIR, "subcatch_tiffs")
EXPORT_SEPARATE_TIFFS = True
WRITE_MASKS_AS = "mask"     # "mask" (0/1) or "id" (id elsewhere nodata)
MAX_EXPORT = None           # set e.g. 200 for testing
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


def tile_dim(dim, preferred=512):
    """Valid GeoTIFF tile size: multiple of 16 and <= dim, else None."""
    if dim < 16:
        return None
    t = min(preferred, dim)
    t = (t // 16) * 16
    return max(16, t)


# ESRI D8 pointer codes:
# 1=E, 2=SE, 4=S, 8=SW, 16=W, 32=NW, 64=N, 128=NE
# For neighbor at (ndy, ndx) relative to current, to flow into current,
# neighbor must point (-ndy, -ndx):
STEP_TO_CODE = {
    ( 0,  1): 1,    # E
    ( 1,  1): 2,    # SE
    ( 1,  0): 4,    # S
    ( 1, -1): 8,    # SW
    ( 0, -1): 16,   # W
    (-1, -1): 32,   # NW
    (-1,  0): 64,   # N
    (-1,  1): 128,  # NE
}


def expand_window(win: Window, pad: int, h: int, w: int) -> Window:
    r0 = max(0, win.row_off - pad)
    c0 = max(0, win.col_off - pad)
    r1 = min(h, win.row_off + win.height + pad)
    c1 = min(w, win.col_off + win.width + pad)
    return Window(c0, r0, c1 - c0, r1 - r0)  # note: Window(col_off, row_off,...)


def compute_pour_raster_blockwise(pntr_tif: str, uaa_km2_tif: str, out_pour_tif: str, thresh_km2: float):
    """
    Build pour-point raster (int32) in one pass:
      - stream = (UAA >= thresh)
      - inflow count computed from ESRI pointer + stream in a 1-cell halo
      - pour points = stream & (inflow==0 or inflow>=2)
      - assign sequential IDs to pour pixels

    Writes out_pour_tif with nodata=0 and values 1..N.
    """
    print("Building pour-point raster (blockwise)...")

    with rasterio.open(pntr_tif) as pntr_ds, rasterio.open(uaa_km2_tif) as uaa_ds:
        if (pntr_ds.width != uaa_ds.width) or (pntr_ds.height != uaa_ds.height) or (pntr_ds.transform != uaa_ds.transform):
            raise RuntimeError("Pointer and UAA grids are not aligned (shape/transform mismatch).")

        h, w = pntr_ds.height, pntr_ds.width
        pntr_nd = pntr_ds.nodata
        uaa_nd = uaa_ds.nodata

        prof = pntr_ds.profile.copy()
        prof.update(
            driver="GTiff",
            dtype="int32",
            count=1,
            nodata=0,
            compress="lzw",
            BIGTIFF="IF_SAFER"
        )

        # safe tiling for large rasters
        bx = tile_dim(w)
        by = tile_dim(h)
        if bx is not None and by is not None:
            prof.update(tiled=True, blockxsize=bx, blockysize=by)
        else:
            prof.update(tiled=False)
            prof.pop("blockxsize", None)
            prof.pop("blockysize", None)

        # write temp then rename (prevents half-written TIFF issues)
        out_path = Path(out_pour_tif)
        tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")
        if tmp_path.exists():
            tmp_path.unlink()
        if out_path.exists():
            out_path.unlink()

        id_counter = 0

        with rasterio.open(tmp_path, "w", **prof) as out_ds:
            for _, win in pntr_ds.block_windows(1):
                win = Window(win.col_off, win.row_off, win.width, win.height)
                big = expand_window(win, pad=1, h=h, w=w)

                pntr_big = pntr_ds.read(1, window=big)
                uaa_big = uaa_ds.read(1, window=big)

                # stream in big window
                stream_big = np.isfinite(uaa_big)
                if uaa_nd is not None:
                    stream_big &= (uaa_big != uaa_nd)
                stream_big &= (uaa_big >= thresh_km2)

                # pointer valid
                pntr_valid_big = np.isfinite(pntr_big)
                if pntr_nd is not None:
                    pntr_valid_big &= (pntr_big != pntr_nd)
                pntr_valid_big &= (pntr_big != 0)

                inflow_big = np.zeros(pntr_big.shape, dtype=np.uint8)

                # compute inflow counts
                for ndy in (-1, 0, 1):
                    for ndx in (-1, 0, 1):
                        if ndy == 0 and ndx == 0:
                            continue
                        needed_code = STEP_TO_CODE.get((-ndy, -ndx))
                        if needed_code is None:
                            continue

                        H, W = pntr_big.shape

                        # dst indices where neighbor exists
                        r0 = max(0, -ndy)
                        r1 = H - max(0, ndy)
                        c0 = max(0, -ndx)
                        c1 = W - max(0, ndx)

                        dst = (slice(r0, r1), slice(c0, c1))
                        src = (slice(r0 + ndy, r1 + ndy), slice(c0 + ndx, c1 + ndx))

                        nbr_ok = (
                            stream_big[src]
                            & pntr_valid_big[src]
                            & (pntr_big[src] == needed_code)
                        )
                        inflow_big[dst] += nbr_ok.astype(np.uint8)

                # crop back to interior window (remove halo)
                r_off = int(win.row_off - big.row_off)
                c_off = int(win.col_off - big.col_off)
                rr = slice(r_off, r_off + int(win.height))
                cc = slice(c_off, c_off + int(win.width))

                stream = stream_big[rr, cc]
                inflow = inflow_big[rr, cc]

                pour_mask = stream & ((inflow == 0) | (inflow >= 2))

                block_out = np.zeros((int(win.height), int(win.width)), dtype=np.int32)
                n = int(pour_mask.sum())
                if n > 0:
                    ids = np.arange(id_counter + 1, id_counter + n + 1, dtype=np.int32)
                    block_out[pour_mask] = ids
                    id_counter += n

                out_ds.write(block_out, 1, window=win)

        tmp_path.replace(out_path)
        print(f"Pour raster written: {out_pour_tif}")
        print(f"Pour points count: {id_counter}")


def export_catchments_cropped_streaming(subcatch_tif: str, export_dir: str, write_as="mask", max_export=None):
    """
    Streaming export:
      Pass 1: scan blocks -> bounding boxes per catchment ID
      Pass 2: write each catchment cropped to bbox as GeoTIFF
    """
    Path(export_dir).mkdir(parents=True, exist_ok=True)

    with rasterio.open(subcatch_tif) as ds:
        nd = ds.nodata
        h, w = ds.height, ds.width
        tfm0 = ds.transform
        prof0 = ds.profile.copy()

        # Pass 1: bounding boxes
        boxes = {}  # id -> [minr, maxr, minc, maxc]
        for _, win in ds.block_windows(1):
            arr = ds.read(1, window=win)

            valid = np.isfinite(arr)
            if nd is not None:
                valid &= (arr != nd)
            valid &= (arr != 0)

            if not valid.any():
                continue

            rr, cc = np.nonzero(valid)
            vals = arr[rr, cc].astype(np.int64)

            r_global = rr + int(win.row_off)
            c_global = cc + int(win.col_off)

            for cid in np.unique(vals):
                m = (vals == cid)
                rmin = int(r_global[m].min()); rmax = int(r_global[m].max())
                cmin = int(c_global[m].min()); cmax = int(c_global[m].max())
                b = boxes.get(int(cid))
                if b is None:
                    boxes[int(cid)] = [rmin, rmax, cmin, cmax]
                else:
                    b[0] = min(b[0], rmin); b[1] = max(b[1], rmax)
                    b[2] = min(b[2], cmin); b[3] = max(b[3], cmax)

        ids = sorted(boxes.keys())
        if max_export is not None and len(ids) > max_export:
            ids = ids[:max_export]

        print(f"Exporting {len(ids)} catchments to: {export_dir}")

        # Output profile for each export
        if write_as.lower() == "mask":
            out_dtype = "uint8"
            out_nodata = 0
        elif write_as.lower() == "id":
            out_dtype = prof0["dtype"]
            out_nodata = nd if nd is not None else 0
        else:
            raise ValueError("write_as must be 'mask' or 'id'")

        for i, cid in enumerate(ids, start=1):
            rmin, rmax, cmin, cmax = boxes[cid]
            r0, r1 = rmin, rmax + 1
            c0, c1 = cmin, cmax + 1

            win = Window(c0, r0, c1 - c0, r1 - r0)
            arr = ds.read(1, window=win)

            if write_as.lower() == "mask":
                out_arr = (arr == cid).astype(np.uint8)
            else:
                out_arr = np.where(arr == cid, arr, out_nodata)

            out_tfm = rasterio.windows.transform(win, tfm0)
            hh, ww = out_arr.shape

            bx = tile_dim(ww)
            by = tile_dim(hh)

            out_prof = prof0.copy()
            out_prof.update(
                driver="GTiff",
                height=hh,
                width=ww,
                transform=out_tfm,
                dtype=out_dtype,
                nodata=out_nodata,
                count=1,
                compress="lzw",
                BIGTIFF="IF_SAFER",
            )
            if bx is not None and by is not None:
                out_prof.update(tiled=True, blockxsize=bx, blockysize=by)
            else:
                out_prof.update(tiled=False)
                out_prof.pop("blockxsize", None)
                out_prof.pop("blockysize", None)

            out_path = Path(export_dir) / f"subcatch_{cid:08d}.tif"
            tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")

            if tmp_path.exists():
                tmp_path.unlink()
            if out_path.exists():
                out_path.unlink()

            with rasterio.open(tmp_path, "w", **out_prof) as out_ds:
                out_ds.write(out_arr.astype(out_dtype), 1)

            tmp_path.replace(out_path)

            if i % 250 == 0:
                print(f"  ... {i}/{len(ids)} written")

        print("Done exporting catchments.")


# -------------------------
# MAIN
# -------------------------
print("0) Choose DEM input...")
dem_in = DEM_90M
if DEM_MASKED:
    dem_in = DEM_MASKED
if DEM_BURNED:
    dem_in = DEM_BURNED
assert os.path.exists(dem_in), f"DEM not found: {dem_in}"
print("DEM:", dem_in)

print("1) Fill depressions (recommended)...")
if not os.path.exists(DEM_FILLED):
    run_wbt(
        [
            "--run=FillDepressions",
            f"--dem={dem_in}",
            f"--output={DEM_FILLED}",
        ],
        wd=OUT_DIR
    )
print("   Wrote:", DEM_FILLED)

print("2) D8 pointer (ESRI) from DEM...")
if not os.path.exists(D8_PNTR):
    run_wbt(
        [
            "--run=D8Pointer",
            f"--dem={DEM_FILLED}",
            f"--output={D8_PNTR}",
        ],
        wd=OUT_DIR
    )
print("   Wrote:", D8_PNTR)

print("3) Flow accumulation (cells) on FULL domain...")
if not os.path.exists(UP_CELLS):
    run_wbt(
        [
            "--run=D8FlowAccumulation",
            f"--input={D8_PNTR}",
            f"--output={UP_CELLS}",
            "--out_type=cells",
            "--pntr",
            "--esri_pntr",
        ],
        wd=OUT_DIR
    )
print("   Wrote:", UP_CELLS)

print("4) Convert upstream cells -> upstream area (km²) ...")
if not os.path.exists(UAA_KM2):
    with rasterio.open(D8_PNTR) as ds_p:
        resx, resy = ds_p.res
        # projected CRS assumed (meters); constant cell area
        cell_km2 = (abs(resx) * abs(resy)) / 1e6

    with rasterio.open(UP_CELLS) as ds:
        prof = ds.profile.copy()
        nd = ds.nodata

        prof.update(driver="GTiff", dtype="float32", count=1, compress="lzw", BIGTIFF="IF_SAFER")

        bx = tile_dim(ds.width)
        by = tile_dim(ds.height)
        if bx is not None and by is not None:
            prof.update(tiled=True, blockxsize=bx, blockysize=by)
        else:
            prof.update(tiled=False)
            prof.pop("blockxsize", None)
            prof.pop("blockysize", None)

        out_path = Path(UAA_KM2)
        tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")
        if tmp_path.exists():
            tmp_path.unlink()
        if out_path.exists():
            out_path.unlink()

        with rasterio.open(tmp_path, "w", **prof) as out_ds:
            for _, win in ds.block_windows(1):
                up = ds.read(1, window=win).astype("float64")
                uaa = up * cell_km2
                if nd is not None:
                    uaa[up == nd] = nd
                out_ds.write(uaa.astype("float32"), 1, window=win)

        tmp_path.replace(out_path)

print("   Wrote:", UAA_KM2)

print("5) Build pour points raster (headwaters + confluences) from UAA>=1000km² ...")
if not os.path.exists(POUR_RASTER):
    compute_pour_raster_blockwise(D8_PNTR, UAA_KM2, POUR_RASTER, THRESH_KM2)
print("   Wrote:", POUR_RASTER)

print("6) Watershed (pointer + pour raster) -> subcatch raster ...")
if not os.path.exists(SUBCATCH_RASTER):
    run_wbt(
        [
            "--run=Watershed",
            f"--d8_pntr={D8_PNTR}",
            f"--pour_pts={POUR_RASTER}",
            f"--output={SUBCATCH_RASTER}",
            "--esri_pntr",
        ],
        wd=OUT_DIR
    )
print("   Wrote:", SUBCATCH_RASTER)

print("7) Raster -> polygons (optional; can be huge) ...")
# Uncomment if needed
# if not os.path.exists(SUBCATCH_GPKG):
#     run_wbt(
#         [
#             "--run=RasterToVectorPolygons",
#             f"--input={SUBCATCH_RASTER}",
#             f"--output={SUBCATCH_GPKG}",
#         ],
#         wd=OUT_DIR
#     )
# print("   Wrote:", SUBCATCH_GPKG)

print("8) Export each catchment as cropped GeoTIFFs ...")
if EXPORT_SEPARATE_TIFFS:
    export_catchments_cropped_streaming(
        SUBCATCH_RASTER,
        EXPORT_DIR,
        write_as=WRITE_MASKS_AS,
        max_export=MAX_EXPORT
    )

print("DONE.")
