#!/usr/bin/env python3
"""
EU subcatchments from a 90m DEM (ETRS89 / LAEA Europe, EPSG:3035), large-raster-safe.

Your 90m DEM:
  E:\pythonProject\Europe\PREP_90M\dem_90m_3035.tif
  cellsize: 90 x 90 m
  rows/cols: 48533 x 55548  (very large)

Workflow (WhiteboxTools + optional ArcGIS river burning):
  0) (Optional) burn rivers into DEM using ArcGIS Pro (arcpy)
  1) Fill depressions (Whitebox FillDepressions)
  2) D8 pointer (Whitebox D8Pointer) -> ESRI pointer (powers of two)
  3) D8 flow accumulation (cells) (Whitebox D8FlowAccumulation)
  4) Convert cells -> upstream area (km²), constant cell area (90m*90m)
  5) Stream mask: UAA_km2 >= THRESH_KM2 (default 1000 km²)
  6) Pour points raster: headwaters + confluences on stream network (blockwise)
  7) Watershed (Whitebox Watershed) -> subcatch raster
  8) Export each catchment as cropped GeoTIFFs (streamed; safe tiling)

Run in PyCharm with ArcGIS Pro Python if you want burning:
  "C:\\Program Files\\ArcGIS\\Pro\\bin\\Python\\envs\\arcgispro-py3\\python.exe" script.py
Otherwise set RIVERS=None and it runs without arcpy.

Notes:
- This is EU-scale; exports can be thousands of files.
- Uses temp-write + rename to avoid partially-written TIFF issues.

Fixes:
- avoids writing projected rivers to .gpkg (ArcPy expects a feature class path)
- uses FileGDB for projected rivers
- uses OIDFieldName for PolylineToRaster value_field
- removes EPSG:3035 strict check; only requires projected meters for km² conversion
"""

import os
import time
from pathlib import Path
import numpy as np
import subprocess

import rasterio
from rasterio.windows import Window

# -------------------------
# INPUTS (EDIT)
# -------------------------
OUT_DIR = r"E:\pythonProject\Europe\OUT_DEM90M_WHITEBOX"
Path(OUT_DIR).mkdir(parents=True, exist_ok=True)

DEM_90M = r"E:\pythonProject\Europe\PREP_90M\dem_90m_3035.tif"

LAND_MASK = None  # optional aligned land mask raster (land=1, water=0), processed with ArcGIS if provided

RIVERS = r"E:\pythonProject\spu\river.shp"  # set None to disable burning
BURN_DEPTH_M = 10.0
RIVER_BUFFER_CELLS = 1

WBT_EXE = r"C:\Users\MRE\anaconda4\envs\Subcatchment\Library\bin\whitebox_tools.exe"
assert os.path.exists(WBT_EXE), f"whitebox_tools.exe not found: {WBT_EXE}"

THRESH_KM2 = 25.0  # stream definition threshold

# Outputs / intermediates
DEM_LAND     = os.path.join(OUT_DIR, "dem_land_only.tif")
DEM_BURNED   = os.path.join(OUT_DIR, "dem_burned.tif")
DEM_FILLED   = os.path.join(OUT_DIR, "dem_filled.tif")
D8_PNTR      = os.path.join(OUT_DIR, "d8_pointer_esri.tif")
UP_CELLS     = os.path.join(OUT_DIR, "upstream_cells.tif")
UAA_KM2      = os.path.join(OUT_DIR, "upstream_area_km2.tif")
POUR_RASTER  = os.path.join(OUT_DIR, "pour_pts.tif")
SUBCATCH_RASTER = os.path.join(OUT_DIR, "subcatch.tif")

EXPORT_DIR = os.path.join(OUT_DIR, "subcatch_tiffs")
EXPORT_SEPARATE_TIFFS = True
WRITE_MASKS_AS = "mask"   # "mask" (0/1) or "id"
MAX_EXPORT = None         # e.g. 200 for testing
# -------------------------


# -------------------------
# Timing helpers
# -------------------------
class Timer:
    def __init__(self):
        self.t0 = time.perf_counter()
        self.marks = []

    def mark(self, label: str):
        t = time.perf_counter()
        self.marks.append((label, t - self.t0))

    def summary(self):
        if not self.marks:
            return
        print("\n--- Timing summary ---")
        prev = 0.0
        for lbl, tsec in self.marks:
            print(f"{lbl:<28} {tsec/60:.2f} min  (Δ {(tsec-prev)/60:.2f} min)")
            prev = tsec


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


# ESRI D8 pointer codes
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
    return Window(c0, r0, c1 - c0, r1 - r0)


def ensure_filegdb(gdb_path: str):
    import arcpy
    folder = os.path.dirname(gdb_path)
    name = os.path.basename(gdb_path)
    if not arcpy.Exists(gdb_path):
        arcpy.management.CreateFileGDB(folder, name)


def burn_rivers_into_dem_arcpy(dem_path, rivers_fc, out_dem, burn_depth_m=10.0, buffer_cells=1):
    """
    Burn rivers into DEM using ArcGIS:
      - project rivers to DEM CRS into a FileGDB feature class
      - polyline->raster aligned to DEM
      - optional Expand to widen channel
      - subtract burn depth where river raster is 1
    """
    import arcpy
    from arcpy.sa import Raster, Con, IsNull, Expand

    arcpy.CheckOutExtension("Spatial")
    arcpy.env.overwriteOutput = True
    arcpy.env.snapRaster = dem_path
    arcpy.env.cellSize = dem_path
    arcpy.env.extent = dem_path
    arcpy.env.workspace = OUT_DIR

    scratch_gdb = os.path.join(OUT_DIR, "_scratch.gdb")
    ensure_filegdb(scratch_gdb)

    if not arcpy.Exists(rivers_fc):
        raise FileNotFoundError(f"RIVERS not found: {rivers_fc}")

    dem_sr = arcpy.Describe(dem_path).spatialReference
    riv_sr = arcpy.Describe(rivers_fc).spatialReference

    rivers_use = rivers_fc
    rivers_proj_fc = os.path.join(scratch_gdb, "rivers_proj")

    # Reproject if needed (or just overwrite projected fc)
    if arcpy.Exists(rivers_proj_fc):
        arcpy.management.Delete(rivers_proj_fc)

    if riv_sr and dem_sr and (riv_sr.factoryCode != dem_sr.factoryCode or riv_sr.name != dem_sr.name):
        arcpy.management.Project(rivers_fc, rivers_proj_fc, dem_sr)
        rivers_use = rivers_proj_fc
    else:
        # still copy into gdb to guarantee support
        arcpy.management.CopyFeatures(rivers_fc, rivers_proj_fc)
        rivers_use = rivers_proj_fc

    oid_field = arcpy.Describe(rivers_use).OIDFieldName

    dem = Raster(dem_path)
    rivers_ras = os.path.join(OUT_DIR, "rivers_ras.tif")
    if arcpy.Exists(rivers_ras):
        arcpy.management.Delete(rivers_ras)

    arcpy.conversion.PolylineToRaster(
        in_features=rivers_use,
        value_field=oid_field,
        out_rasterdataset=rivers_ras,
        cell_assignment="MAXIMUM_LENGTH",
        priority_field="NONE",
        cellsize=dem.meanCellWidth
    )

    r = Raster(rivers_ras)
    r01 = Con(IsNull(r), 0, 1)

    if buffer_cells and int(buffer_cells) > 0:
        r01 = Expand(r01, int(buffer_cells), [1])

    burned = dem - (r01 * float(burn_depth_m))

    # Safe overwrite
    if os.path.exists(out_dem):
        try:
            os.remove(out_dem)
        except OSError:
            pass
    burned.save(out_dem)
    return out_dem


def apply_land_mask_arcpy(dem_path, land_mask_path, out_dem_path):
    import arcpy
    from arcpy.sa import Raster, SetNull

    arcpy.CheckOutExtension("Spatial")
    arcpy.env.overwriteOutput = True
    arcpy.env.snapRaster = dem_path
    arcpy.env.cellSize = dem_path
    arcpy.env.extent = dem_path
    arcpy.env.workspace = OUT_DIR

    if not arcpy.Exists(land_mask_path):
        raise FileNotFoundError(f"LAND_MASK not found: {land_mask_path}")

    dem = Raster(dem_path)
    land = Raster(land_mask_path)
    dem_land = SetNull(land == 0, dem)

    if os.path.exists(out_dem_path):
        try:
            os.remove(out_dem_path)
        except OSError:
            pass
    dem_land.save(out_dem_path)
    return out_dem_path


def compute_pour_raster_blockwise(pntr_tif: str, uaa_km2_tif: str, out_pour_tif: str, thresh_km2: float):
    print("Building pour-point raster (blockwise)...")

    with rasterio.open(pntr_tif) as pntr_ds, rasterio.open(uaa_km2_tif) as uaa_ds:
        if (pntr_ds.width != uaa_ds.width) or (pntr_ds.height != uaa_ds.height) or (pntr_ds.transform != uaa_ds.transform):
            raise RuntimeError("Pointer and UAA grids are not aligned (shape/transform mismatch).")

        h, w = pntr_ds.height, pntr_ds.width
        pntr_nd = pntr_ds.nodata
        uaa_nd = uaa_ds.nodata

        prof = pntr_ds.profile.copy()
        prof.update(driver="GTiff", dtype="int32", count=1, nodata=0, compress="lzw", BIGTIFF="IF_SAFER")

        bx = tile_dim(w)
        by = tile_dim(h)
        if bx is not None and by is not None:
            prof.update(tiled=True, blockxsize=bx, blockysize=by)
        else:
            prof.update(tiled=False)
            prof.pop("blockxsize", None)
            prof.pop("blockysize", None)

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

                stream_big = np.isfinite(uaa_big)
                if uaa_nd is not None:
                    stream_big &= (uaa_big != uaa_nd)
                stream_big &= (uaa_big >= thresh_km2)

                pntr_valid_big = np.isfinite(pntr_big)
                if pntr_nd is not None:
                    pntr_valid_big &= (pntr_big != pntr_nd)
                pntr_valid_big &= (pntr_big != 0)

                inflow_big = np.zeros(pntr_big.shape, dtype=np.uint8)

                for ndy in (-1, 0, 1):
                    for ndx in (-1, 0, 1):
                        if ndy == 0 and ndx == 0:
                            continue
                        needed_code = STEP_TO_CODE.get((-ndy, -ndx))
                        if needed_code is None:
                            continue

                        H, W = pntr_big.shape
                        r0 = max(0, -ndy)
                        r1 = H - max(0, ndy)
                        c0 = max(0, -ndx)
                        c1 = W - max(0, ndx)

                        dst = (slice(r0, r1), slice(c0, c1))
                        src = (slice(r0 + ndy, r1 + ndy), slice(c0 + ndx, c1 + ndx))

                        nbr_ok = stream_big[src] & pntr_valid_big[src] & (pntr_big[src] == needed_code)
                        inflow_big[dst] += nbr_ok.astype(np.uint8)

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
    Path(export_dir).mkdir(parents=True, exist_ok=True)

    with rasterio.open(subcatch_tif) as ds:
        nd = ds.nodata
        tfm0 = ds.transform
        prof0 = ds.profile.copy()

        boxes = {}
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
            win = Window(cmin, rmin, (cmax - cmin + 1), (rmax - rmin + 1))
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
T = Timer()

print("0) Checking DEM...")
assert os.path.exists(DEM_90M), f"DEM not found: {DEM_90M}"

with rasterio.open(DEM_90M) as ds:
    print("DEM:", DEM_90M)
    print("CRS:", ds.crs)
    print("RES:", ds.res)
    print("Rows/Cols:", ds.height, ds.width)
    if ds.crs is None or not ds.crs.is_projected:
        raise RuntimeError("DEM CRS must be projected (meters) for constant km² cell-area conversion.")
T.mark("DEM check")

dem_in = DEM_90M

if LAND_MASK:
    print("0b) Applying LAND_MASK (ArcGIS)...")
    dem_in = apply_land_mask_arcpy(dem_in, LAND_MASK, DEM_LAND)
    print("Masked DEM:", dem_in)
    T.mark("Apply land mask")

if RIVERS:
    print("0c) Burning rivers into DEM (ArcGIS)...")
    dem_in = burn_rivers_into_dem_arcpy(
        dem_in, RIVERS, DEM_BURNED,
        burn_depth_m=BURN_DEPTH_M,
        buffer_cells=RIVER_BUFFER_CELLS
    )
    print("Burned DEM:", dem_in)
    T.mark("Burn rivers")

print("1) Fill depressions...")
if not os.path.exists(DEM_FILLED):
    run_wbt(["--run=FillDepressions", f"--dem={dem_in}", f"--output={DEM_FILLED}"], wd=OUT_DIR)
print("   Wrote:", DEM_FILLED)
T.mark("Fill depressions")

print("2) D8 pointer (ESRI) from DEM...")
if not os.path.exists(D8_PNTR):
    run_wbt(["--run=D8Pointer", f"--dem={DEM_FILLED}", f"--output={D8_PNTR}"], wd=OUT_DIR)
print("   Wrote:", D8_PNTR)
T.mark("D8 pointer")

print("3) Flow accumulation (cells) on full domain...")
if not os.path.exists(UP_CELLS):
    run_wbt(
        ["--run=D8FlowAccumulation", f"--input={D8_PNTR}", f"--output={UP_CELLS}",
         "--out_type=cells", "--pntr", "--esri_pntr"],
        wd=OUT_DIR
    )
print("   Wrote:", UP_CELLS)
T.mark("Flow accumulation")

print("4) Convert upstream cells -> upstream area (km²)...")
with rasterio.open(D8_PNTR) as ds:
    resx, resy = ds.res
cell_km2 = (abs(resx) * abs(resy)) / 1e6
print("   cell_km2 =", cell_km2)

if not os.path.exists(UAA_KM2):
    with rasterio.open(UP_CELLS) as ds:
        prof = ds.profile.copy()
        nd = ds.nodata

        prof.update(driver="GTiff", dtype="float32", count=1, nodata=nd, compress="lzw", BIGTIFF="IF_SAFER")

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
T.mark("Cells -> km²")

print(f"5) Build pour points raster (UAA >= {THRESH_KM2} km²)...")
if not os.path.exists(POUR_RASTER):
    compute_pour_raster_blockwise(D8_PNTR, UAA_KM2, POUR_RASTER, THRESH_KM2)
print("   Wrote:", POUR_RASTER)
T.mark("Pour points raster")

print("6) Watershed -> subcatch raster...")
if not os.path.exists(SUBCATCH_RASTER):
    run_wbt(
        ["--run=Watershed", f"--d8_pntr={D8_PNTR}", f"--pour_pts={POUR_RASTER}", f"--output={SUBCATCH_RASTER}", "--esri_pntr"],
        wd=OUT_DIR
    )
print("   Wrote:", SUBCATCH_RASTER)
T.mark("Watershed")

print("7) Export catchments (cropped GeoTIFFs)...")
if EXPORT_SEPARATE_TIFFS:
    export_catchments_cropped_streaming(SUBCATCH_RASTER, EXPORT_DIR, write_as=WRITE_MASKS_AS, max_export=MAX_EXPORT)
T.mark("Export catchments")

print("DONE.")
T.summary()
