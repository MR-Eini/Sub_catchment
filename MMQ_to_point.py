#!/usr/bin/env python3
import os
import arcpy
from arcpy.sa import Raster, SetNull

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True

in_ras  = r"E:\pythonProject\Europe\MMQ_mamad.tif"
out_pts = r"E:\pythonProject\Europe\MMQ_mamad_points.shp"   # or .gpkg/.gdb FC

# OPTIONAL (recommended):
# 1) remove negatives -> NoData
# 2) remove zeros -> NoData (uncomment if you don't want points for 0)
r = Raster(in_ras)
r_clean = SetNull(r < 0, r)
# r_clean = SetNull((r_clean == 0), r_clean)

tmp_ras = r"E:\pythonProject\Europe\_tmp_MMQ_clean.tif"
r_clean.save(tmp_ras)

# Raster cell centers -> points; values go to field "grid_code"
arcpy.conversion.RasterToPoint(tmp_ras, out_pts, "VALUE")

# Create MMQ field and copy values from grid_code
arcpy.management.AddField(out_pts, "MMQ", "DOUBLE")
arcpy.management.CalculateField(out_pts, "MMQ", "!grid_code!", "PYTHON3")

# Optional: drop the default field
arcpy.management.DeleteField(out_pts, ["grid_code"])

print("DONE:", out_pts)
