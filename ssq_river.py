import os
import arcpy
from arcpy.sa import Raster, SetNull, Thin, Int, ZonalStatisticsAsTable

arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")

ssq_raster = r"E:\pythonProject\Europe\SSQ.tif"

out_dir  = r"E:\pythonProject\Europe"
gdb_name = "SSQ_rivers.gdb"
out_gdb  = os.path.join(out_dir, gdb_name)
out_fc   = os.path.join(out_gdb, "ssq_rivers")

# OPTIONAL (recommended): if you have D8 flow direction raster, set it here
flow_dir_raster = None
# flow_dir_raster = r"E:\pythonProject\Europe\flowdir.tif"

os.makedirs(out_dir, exist_ok=True)
if not arcpy.Exists(out_gdb):
    arcpy.management.CreateFileGDB(out_dir, gdb_name)
if arcpy.Exists(out_fc):
    arcpy.management.Delete(out_fc)

scratch = arcpy.env.scratchGDB or out_gdb

r = Raster(ssq_raster)

# Stream mask: 1 where discharge > 0, else NoData
stream_ras = SetNull(r <= 0, 1)

# --- Raster -> Polyline ---
if flow_dir_raster and arcpy.Exists(flow_dir_raster):
    # Best method (uses flow direction)
    arcpy.sa.StreamToFeature(stream_ras, Raster(flow_dir_raster), out_fc, "SIMPLIFY")
else:
    # Fallback: thin then raster-to-polyline
    thin_ras = Thin(stream_ras, "ZERO", "NO_FILTER")
    thin_int = Int(thin_ras)  # RasterToPolyline prefers integer rasters

    # Correct parameter order:
    # RasterToPolyline(in_raster, out_fc, background_value, minimum_dangle_length, simplify)
    arcpy.conversion.RasterToPolyline(thin_int, out_fc, "NODATA", 0, "NO_SIMPLIFY")

# --- Add discharge attributes to polylines (mean/max over intersecting cells) ---
zone_field = "RID"
if zone_field.upper() not in {f.name.upper() for f in arcpy.ListFields(out_fc)}:
    arcpy.management.AddField(out_fc, zone_field, "LONG")
    arcpy.management.CalculateField(out_fc, zone_field, "!OBJECTID!", "PYTHON3")

zs_tbl = arcpy.CreateUniqueName("zs_ssq", scratch)
ZonalStatisticsAsTable(out_fc, zone_field, r, zs_tbl, "DATA", "ALL")

arcpy.management.JoinField(out_fc, zone_field, zs_tbl, zone_field, ["MEAN", "MAX"])

for fld, src in [("Q_MEAN", "MEAN"), ("Q_MAX", "MAX")]:
    if fld.upper() not in {f.name.upper() for f in arcpy.ListFields(out_fc)}:
        arcpy.management.AddField(out_fc, fld, "DOUBLE")
    arcpy.management.CalculateField(out_fc, fld, f"!{src}!", "PYTHON3")

# Optional cleanup of join fields
try:
    arcpy.management.DeleteField(out_fc, ["MEAN", "MAX"])
except Exception:
    pass
try:
    arcpy.management.Delete(zs_tbl)
except Exception:
    pass

print(f"Done. Output polyline: {out_fc}")
print("Fields: Q_MEAN, Q_MAX")
