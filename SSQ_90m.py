import os, arcpy
from arcpy.sa import Raster
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True

BASE = r"E:\pythonProject\Europe\OUT_ARCGIS_90M"
SSQ_IN = r"E:\pythonProject\Europe\SSQ.tif"

SNAP = os.path.join(BASE, "streamlink.tif")   # 90m template
OUT  = r"E:\pythonProject\Europe\DELTAQ_90M"
os.makedirs(OUT, exist_ok=True)

SSQ_3035 = os.path.join(OUT, "SSQ_90m_3035.tif")

arcpy.env.snapRaster = SNAP
arcpy.env.cellSize   = SNAP
arcpy.env.extent     = SNAP

# If SSQ is continuous discharge/runoff: use BILINEAR.
# If SSQ is categorical: use NEAREST.
arcpy.management.ProjectRaster(
    in_raster=SSQ_IN,
    out_raster=SSQ_3035,
    out_coor_system=arcpy.Describe(SNAP).spatialReference,
    resampling_type="BILINEAR",
    cell_size=arcpy.Describe(SNAP).meanCellWidth
)

print("Wrote aligned SSQ:", SSQ_3035)
