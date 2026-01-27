#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import csv
import arcpy
from arcpy.sa import Raster, SetNull, ZonalStatisticsAsTable

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True

# -----------------------------
# INPUTS
# -----------------------------
RATIO_IN = r"Y:\download\Mamad\ratio_p5_p95.tif"
SUBCATCH = r"E:\pythonProject\Europe\OUT_ARCGIS_25M\subcatch_streamlink.tif"  # zone raster (Value=ID)

OUT_DIR  = r"E:\pythonProject\Europe\RATIO_P5_P95_SUBCATCH"
os.makedirs(OUT_DIR, exist_ok=True)

# -----------------------------
# SETTINGS
# -----------------------------
# If input CRS is unknown, we will DefineProjection using this EPSG
# (Heuristic will override to 4326 if raster looks like degrees, or 3035 if looks like meters)
ASSUMED_EPSG_IF_UNKNOWN = 4326

# Treat very small values as NoData-like
# Start with 1e-6; if too aggressive use 1e-7; if not removing use 1e-5
NODATA_EPS = 1e-6

MAKE_SHP = True  # optional polygons + dissolve + join + export shapefile

# -----------------------------
# WORKSPACE
# -----------------------------
GDB = os.path.join(OUT_DIR, "ratio_subcatch.gdb")
if not arcpy.Exists(GDB):
    arcpy.management.CreateFileGDB(OUT_DIR, "ratio_subcatch.gdb")
arcpy.env.workspace = GDB
arcpy.env.scratchWorkspace = GDB

# alignment template
arcpy.env.snapRaster = SUBCATCH
arcpy.env.cellSize   = SUBCATCH
arcpy.env.extent     = SUBCATCH

# -----------------------------
# OUTPUTS
# -----------------------------
ALN_DIR        = os.path.join(OUT_DIR, "aligned")
os.makedirs(ALN_DIR, exist_ok=True)

RATIO_PROJ     = os.path.join(ALN_DIR, "ratio_proj_3035.tif")
RATIO_ALIGNED  = os.path.join(ALN_DIR, "ratio_aligned_3035_90m.tif")

RATIO_CLEAN    = os.path.join(GDB, "ratio_clean")

ZONAL_TBL      = os.path.join(GDB, "ratio_mean_tbl")         # Value + MEAN
ALL_IDS_TBL    = os.path.join(GDB, "all_subcatch_ids")       # from SUBCATCH RAT
FINAL_TBL      = os.path.join(GDB, "ratio_mean_by_subcatch") # ALL IDs + RATIO_MEAN

CSV_OUT        = os.path.join(OUT_DIR, "ratio_mean_by_subcatch.csv")

# optional polygons
POLY_RAW       = os.path.join(GDB, "subcatch_poly_raw")
POLY_DISS      = os.path.join(GDB, "subcatch_poly_diss")
POLY_JOIN      = os.path.join(GDB, "subcatch_poly_join")
SHP_OUT        = os.path.join(OUT_DIR, "subcatch_with_ratio_mean.shp")


def is_unknown_sr(sr):
    try:
        if sr is None:
            return True
        # common cases for "unknown"
        return (sr.name is None) or (sr.name.lower() in ("unknown", "unknown coordinate system")) or (getattr(sr, "factoryCode", 0) in (0, None))
    except Exception:
        return True


def infer_epsg_from_extent(ras_path):
    """Heuristic: if extent looks like degrees -> 4326, if looks like meters in Europe -> 3035."""
    d = arcpy.Describe(ras_path)
    e = d.extent
    # degrees-ish
    if (-180 <= e.XMin <= 180 and -180 <= e.XMax <= 180 and -90 <= e.YMin <= 90 and -90 <= e.YMax <= 90):
        return 4326
    # LAEA meters-ish (EU extent typically millions)
    if (abs(e.XMax) > 100000 and abs(e.YMax) > 100000):
        return 3035
    return ASSUMED_EPSG_IF_UNKNOWN


def define_projection_if_needed(ras_path):
    d = arcpy.Describe(ras_path)
    sr = d.spatialReference
    if not is_unknown_sr(sr):
        return

    epsg = infer_epsg_from_extent(ras_path)
    print(f"Input CRS is undefined. Defining projection as EPSG:{epsg} (metadata only)...")
    arcpy.management.DefineProjection(ras_path, arcpy.SpatialReference(epsg))


def safe_csv_path(path):
    if not os.path.exists(path):
        return path
    root, ext = os.path.splitext(path)
    for i in range(1, 1000):
        cand = f"{root}_{i}{ext}"
        if not os.path.exists(cand):
            return cand
    raise RuntimeError("No free CSV filename found.")


def export_table_to_csv(table, out_csv, fields):
    out_csv = safe_csv_path(out_csv)
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(fields)
        with arcpy.da.SearchCursor(table, fields) as cur:
            for row in cur:
                w.writerow(row)
    print("CSV:", out_csv)


def main():
    if not os.path.exists(RATIO_IN):
        raise FileNotFoundError(RATIO_IN)

    print("RATIO_IN:", RATIO_IN)
    print("SUBCATCH:", SUBCATCH)

    # 0) Ensure ratio has a defined CRS
    define_projection_if_needed(RATIO_IN)

    # Template info
    dtemp = arcpy.Describe(SUBCATCH)
    out_sr = dtemp.spatialReference
    cellx = abs(float(dtemp.meanCellWidth))
    celly = abs(float(dtemp.meanCellHeight))
    cellsize = f"{cellx} {celly}"

    print("Template CRS:", out_sr.name)
    print("Template cell:", cellx, celly)

    # 1) Project to template CRS
    print("\n1) Projecting ratio to template CRS...")
    arcpy.env.extent = dtemp.extent
    arcpy.management.ProjectRaster(
        in_raster=RATIO_IN,
        out_raster=RATIO_PROJ,
        out_coor_system=out_sr,
        resampling_type="BILINEAR",
        cell_size=cellsize
    )

    # 2) Snap/resample to exact grid origin
    print("2) Resampling with snap to match template grid origin...")
    arcpy.env.snapRaster = SUBCATCH
    arcpy.env.cellSize   = SUBCATCH
    arcpy.env.extent     = SUBCATCH
    arcpy.management.Resample(
        in_raster=RATIO_PROJ,
        out_raster=RATIO_ALIGNED,
        cell_size=cellsize,
        resampling_type="BILINEAR"
    )
    print("Aligned ratio:", RATIO_ALIGNED)

    # 3) Clean near-zero NoData-like values
    print(f"\n3) Cleaning ratio (<= {NODATA_EPS} -> NoData)...")
    r = Raster(RATIO_ALIGNED)
    r_clean = SetNull(r <= NODATA_EPS, r)
    r_clean.save(RATIO_CLEAN)

    # 4) Zonal mean (may omit zones that are entirely NoData)
    print("\n4) ZonalStatisticsAsTable (MEAN) by subcatch ID...")
    if arcpy.Exists(ZONAL_TBL): arcpy.management.Delete(ZONAL_TBL)
    ZonalStatisticsAsTable(
        in_zone_data=SUBCATCH,
        zone_field="Value",
        in_value_raster=RATIO_CLEAN,
        out_table=ZONAL_TBL,
        ignore_nodata="DATA",
        statistics_type="MEAN"
    )

    # Rename MEAN -> RATIO_MEAN
    mf = [f.name for f in arcpy.ListFields(ZONAL_TBL) if f.name.upper() == "MEAN"]
    if mf:
        arcpy.management.AlterField(ZONAL_TBL, mf[0], "RATIO_MEAN", "RATIO_MEAN")

    # 5) Ensure ALL IDs present: build RAT from SUBCATCH and left-join
    print("\n5) Building master ID table from SUBCATCH RAT and joining mean...")
    arcpy.management.BuildRasterAttributeTable(SUBCATCH, "Overwrite")

    if arcpy.Exists(ALL_IDS_TBL): arcpy.management.Delete(ALL_IDS_TBL)
    arcpy.management.CopyRows(SUBCATCH, ALL_IDS_TBL)  # copies RAT

    if "Value" not in [f.name for f in arcpy.ListFields(ALL_IDS_TBL)]:
        raise RuntimeError("Could not find 'Value' in SUBCATCH raster attribute table.")
    arcpy.management.AlterField(ALL_IDS_TBL, "Value", "SUBC_ID", "SUBC_ID")

    if arcpy.Exists(FINAL_TBL): arcpy.management.Delete(FINAL_TBL)
    arcpy.management.CopyRows(ALL_IDS_TBL, FINAL_TBL)
    arcpy.management.JoinField(FINAL_TBL, "SUBC_ID", ZONAL_TBL, "Value", ["RATIO_MEAN"])

    export_table_to_csv(FINAL_TBL, CSV_OUT, ["SUBC_ID", "RATIO_MEAN"])
    print("Final table:", FINAL_TBL)

    # 6) Optional polygons + dissolve (fix duplicates) + join + shapefile
    if MAKE_SHP:
        print("\n6) Optional: RasterToPolygon + Dissolve + Join + Export shapefile...")
        for fc in [POLY_RAW, POLY_DISS, POLY_JOIN]:
            if arcpy.Exists(fc):
                arcpy.management.Delete(fc)

        arcpy.conversion.RasterToPolygon(SUBCATCH, POLY_RAW, "NO_SIMPLIFY", "VALUE")
        gc = [f.name for f in arcpy.ListFields(POLY_RAW) if f.name.lower() == "gridcode"]
        if gc:
            arcpy.management.AlterField(POLY_RAW, gc[0], "SUBC_ID", "SUBC_ID")

        arcpy.management.Dissolve(POLY_RAW, POLY_DISS, ["SUBC_ID"], multi_part="MULTI_PART")
        arcpy.management.CopyFeatures(POLY_DISS, POLY_JOIN)
        arcpy.management.JoinField(POLY_JOIN, "SUBC_ID", FINAL_TBL, "SUBC_ID", ["RATIO_MEAN"])

        if arcpy.Exists(SHP_OUT):
            arcpy.management.Delete(SHP_OUT)
        arcpy.conversion.FeatureClassToFeatureClass(POLY_JOIN, OUT_DIR, os.path.basename(SHP_OUT))
        print("SHP:", SHP_OUT)

    print("\nDONE.")


if __name__ == "__main__":
    main()
