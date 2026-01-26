#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Export final outputs from existing GDB (no recomputation).
- Exports mmq_flux_by_slink table to CSV (handles locked existing CSV by writing a new name)
- Optional: creates dissolved subcatchment polygons and joins flux, exports shapefile
- Optional: creates flux raster from polygons (DOUBLE) via PolygonToRaster
"""

import os
import arcpy
import csv

arcpy.env.overwriteOutput = True

# -----------------------------
# INPUTS
# -----------------------------
GDB       = r"E:\pythonProject\Europe\MMQ_FLUX_25M\mmq_flux.gdb"
FLUX_TBL  = os.path.join(GDB, "mmq_flux_by_slink")

SUBCATCH  = r"E:\pythonProject\Europe\OUT_ARCGIS_25M\subcatch_streamlink.tif"
SLINK     = r"E:\pythonProject\Europe\OUT_ARCGIS_25M\streamlink.tif"

OUT_DIR   = r"E:\pythonProject\Europe\MMQ_FLUX_25M\EXPORTS"
os.makedirs(OUT_DIR, exist_ok=True)

# What to create
EXPORT_CSV       = True
EXPORT_SHP       = True
EXPORT_FLUX_RAS  = True

# Prefer these fields if present
PREFERRED_FIELDS = ["SLINK_ID", "MMQ_in", "MMQ_out", "MMQ_in_fix", "MMQ_out_fix", "MMQ_flux", "FLAG"]

# -----------------------------
# OUTPUTS
# -----------------------------
CSV_BASE   = os.path.join(OUT_DIR, "mmq_flux_by_slink.csv")
POLY_RAW   = os.path.join(GDB, "subcatch_poly_raw_export")
POLY_DISS  = os.path.join(GDB, "subcatch_poly_diss_export")
POLY_JOIN  = os.path.join(GDB, "subcatch_poly_join_export")
SHP_OUT    = os.path.join(OUT_DIR, "subcatch_with_mmq_flux.shp")
FLUX_RAS   = os.path.join(OUT_DIR, "MMQ_flux_subcatch.tif")


def safe_csv_path(path):
    """If path exists and is locked, generate a new filename."""
    if not os.path.exists(path):
        return path
    root, ext = os.path.splitext(path)
    for i in range(1, 1000):
        candidate = f"{root}_{i}{ext}"
        if not os.path.exists(candidate):
            return candidate
    raise RuntimeError("Could not find a free CSV filename.")


def export_table_to_csv(table, out_csv, fields=None):
    """Write CSV using python (no ArcGIS delete/lock issues)."""
    all_fields = [f.name for f in arcpy.ListFields(table) if f.type not in ("Geometry", "Blob", "Raster")]
    if fields is None:
        use_fields = all_fields
    else:
        use_fields = [f for f in fields if f in all_fields]
        if not use_fields:
            use_fields = all_fields

    out_csv = safe_csv_path(out_csv)

    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(use_fields)
        with arcpy.da.SearchCursor(table, use_fields) as cur:
            for row in cur:
                w.writerow(row)

    print("CSV written:", out_csv)
    print("Fields:", use_fields)


def rename_field_case_insensitive(fc, old_lower, new_name):
    fields = arcpy.ListFields(fc)
    match = [f.name for f in fields if f.name.lower() == old_lower.lower()]
    if match:
        arcpy.management.AlterField(fc, match[0], new_name, new_name)


def make_joined_polygons_and_export():
    # Align outputs to SLINK grid
    arcpy.env.snapRaster = SLINK
    arcpy.env.cellSize   = SLINK
    arcpy.env.extent     = SLINK

    # 1) Raster -> Polygon
    if arcpy.Exists(POLY_RAW): arcpy.management.Delete(POLY_RAW)
    arcpy.conversion.RasterToPolygon(SUBCATCH, POLY_RAW, "NO_SIMPLIFY", "VALUE")
    rename_field_case_insensitive(POLY_RAW, "gridcode", "SLINK_ID")

    # 2) Dissolve to ensure ONE (multipart) feature per SLINK_ID
    if arcpy.Exists(POLY_DISS): arcpy.management.Delete(POLY_DISS)
    arcpy.management.Dissolve(POLY_RAW, POLY_DISS, ["SLINK_ID"], multi_part="MULTI_PART")

    # 3) Copy + Join flux
    if arcpy.Exists(POLY_JOIN): arcpy.management.Delete(POLY_JOIN)
    arcpy.management.CopyFeatures(POLY_DISS, POLY_JOIN)

    # join only fields that exist
    tbl_fields = [f.name for f in arcpy.ListFields(FLUX_TBL)]
    join_fields = [f for f in PREFERRED_FIELDS if f in tbl_fields and f != "SLINK_ID"]
    if "MMQ_flux" not in join_fields and "MMQ_flux" in tbl_fields:
        join_fields.append("MMQ_flux")

    arcpy.management.JoinField(POLY_JOIN, "SLINK_ID", FLUX_TBL, "SLINK_ID", join_fields)
    print("Joined polygons in GDB:", POLY_JOIN)

    # 4) Export shapefile
    if arcpy.Exists(SHP_OUT): arcpy.management.Delete(SHP_OUT)
    arcpy.conversion.FeatureClassToFeatureClass(POLY_JOIN, OUT_DIR, os.path.basename(SHP_OUT))
    print("Shapefile written:", SHP_OUT)

    return POLY_JOIN


def polygon_to_flux_raster(poly_fc):
    # Use SLINK template for exact alignment
    arcpy.env.snapRaster = SLINK
    arcpy.env.cellSize   = SLINK
    arcpy.env.extent     = SLINK

    # Choose field
    fields = [f.name for f in arcpy.ListFields(poly_fc)]
    value_field = "MMQ_flux"
    if value_field not in fields:
        # fallback to *_fix if present
        if "MMQ_flux_fix" in fields:
            value_field = "MMQ_flux_fix"
        else:
            raise RuntimeError("No MMQ_flux or MMQ_flux_fix field found in polygon FC.")

    # Delete existing output if possible
    if os.path.exists(FLUX_RAS):
        try:
            os.remove(FLUX_RAS)
        except Exception:
            pass

    arcpy.conversion.PolygonToRaster(
        in_features=poly_fc,
        value_field=value_field,
        out_rasterdataset=FLUX_RAS,
        cell_assignment="CELL_CENTER",
        priority_field="NONE",
        cellsize=SLINK
    )
    print("Flux raster written:", FLUX_RAS, f"(field={value_field})")


def main():
    if not arcpy.Exists(GDB):
        raise FileNotFoundError(GDB)
    if not arcpy.Exists(FLUX_TBL):
        raise FileNotFoundError(FLUX_TBL)

    arcpy.env.workspace = GDB
    arcpy.env.scratchWorkspace = GDB

    if EXPORT_CSV:
        export_table_to_csv(FLUX_TBL, CSV_BASE, fields=PREFERRED_FIELDS)

    poly_fc = None
    if EXPORT_SHP or EXPORT_FLUX_RAS:
        poly_fc = make_joined_polygons_and_export()

    if EXPORT_FLUX_RAS:
        polygon_to_flux_raster(poly_fc)

    print("DONE.")


if __name__ == "__main__":
    main()
