#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import arcpy

arcpy.env.overwriteOutput = True

# -----------------------------
# INPUTS
# -----------------------------
SUBCATCH_RAS = r"E:\pythonProject\Europe\OUT_ARCGIS_25M\subcatch_streamlink.tif"
CSV_FLUX     = r"E:\pythonProject\Europe\MMQ_FLUX_25M\mmq_flux_by_slink.csv"  # has SLINK_ID, MMQ_in, MMQ_out, MMQ_flux

OUT_DIR      = r"E:\pythonProject\Europe\MMQ_FLUX_25M\POLY_OUT"
os.makedirs(OUT_DIR, exist_ok=True)

# Final shapefile
OUT_SHP      = os.path.join(OUT_DIR, "subcatch_with_mmq_flux.shp")

# Use a FileGDB for intermediate tables (recommended)
GDB = os.path.join(OUT_DIR, "work.gdb")
if not arcpy.Exists(GDB):
    arcpy.management.CreateFileGDB(OUT_DIR, "work.gdb")

POLY_FC = os.path.join(GDB, "subcatch_poly")
TBL_FC  = os.path.join(GDB, "mmq_flux_tbl")

def field_names(tbl):
    return [f.name for f in arcpy.ListFields(tbl)]

def main():
    # 1) Raster -> Polygon (subcatch IDs)
    #    ArcGIS writes the ID field as "gridcode" in the output polygons.
    if arcpy.Exists(POLY_FC):
        arcpy.management.Delete(POLY_FC)

    arcpy.conversion.RasterToPolygon(
        in_raster=SUBCATCH_RAS,
        out_polygon_features=POLY_FC,
        simplify="NO_SIMPLIFY",
        raster_field="VALUE"
    )
    print("Created polygons:", POLY_FC)

    # 2) CSV -> Table in GDB (more reliable than joining directly to CSV)
    if arcpy.Exists(TBL_FC):
        arcpy.management.Delete(TBL_FC)

    arcpy.conversion.TableToTable(CSV_FLUX, GDB, os.path.basename(TBL_FC))
    print("Imported CSV table:", TBL_FC)

    # 3) Ensure join field types match
    #    Polygon join field = "gridcode" (LONG). CSV key expected = "SLINK_ID".
    join_key = "SLINK_ID"
    if join_key not in field_names(TBL_FC):
        raise RuntimeError(f"CSV table does not contain '{join_key}'. Fields: {field_names(TBL_FC)}")

    f = [x for x in arcpy.ListFields(TBL_FC) if x.name == join_key][0]
    if f.type.upper() == "STRING":
        # create numeric key for join
        if "SLINK_ID_L" not in field_names(TBL_FC):
            arcpy.management.AddField(TBL_FC, "SLINK_ID_L", "LONG")
        arcpy.management.CalculateField(TBL_FC, "SLINK_ID_L", "int(!SLINK_ID!)", "PYTHON3")
        tbl_join_field = "SLINK_ID_L"
    else:
        tbl_join_field = "SLINK_ID"

    # 4) Join flux columns into polygons (permanent join)
    #    Keep only needed fields
    keep_fields = []
    for nm in ["MMQ_in", "MMQ_out", "MMQ_flux"]:
        if nm in field_names(TBL_FC):
            keep_fields.append(nm)

    if not keep_fields:
        raise RuntimeError(f"Table has no expected flux fields. Fields: {field_names(TBL_FC)}")

    arcpy.management.JoinField(
        in_data=POLY_FC,
        in_field="gridcode",
        join_table=TBL_FC,
        join_field=tbl_join_field,
        fields=keep_fields
    )
    print("Joined fields:", keep_fields)

    # 5) Export to shapefile
    if arcpy.Exists(OUT_SHP):
        arcpy.management.Delete(OUT_SHP)

    arcpy.conversion.FeatureClassToFeatureClass(
        in_features=POLY_FC,
        out_path=OUT_DIR,
        out_name=os.path.basename(OUT_SHP)
    )
    print("DONE shapefile:", OUT_SHP)

if __name__ == "__main__":
    main()
