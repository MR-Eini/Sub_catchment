#!/usr/bin/env python3
import os
import arcpy
from arcpy.sa import Raster, SetNull

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True

MMQ_IN   = r"E:\pythonProject\Europe\MMQ_mamad.tif"  # WGS84
TEMPLATE = r"E:\pythonProject\Europe\OUT_ARCGIS_25M\streamlink.tif"  # EPSG:3035 grid

OUT_DIR = r"E:\pythonProject\Europe\MMQ_ALIGNED_25M"
os.makedirs(OUT_DIR, exist_ok=True)

MMQ_CLEAN_GCS = os.path.join(OUT_DIR, "MMQ_mamad_clean_wgs84.tif")
MMQ_PROJ_RAW  = os.path.join(OUT_DIR, "MMQ_mamad_3035_raw.tif")
MMQ_ALIGNED   = os.path.join(OUT_DIR, "MMQ_mamad_3035_25m.tif")
MMQ_FINAL     = os.path.join(OUT_DIR, "MMQ_mamad_3035_25m_clean.tif")

# 1) Clean in original CRS FIRST (avoid bilinear smearing of -9999)
mmq = Raster(MMQ_IN)
mmq_clean = SetNull((mmq == -9999) | (mmq < 0), mmq)
mmq_clean.save(MMQ_CLEAN_GCS)

# 2) Project to the template CRS
in_sr  = arcpy.Describe(MMQ_CLEAN_GCS).spatialReference
out_sr = arcpy.Describe(TEMPLATE).spatialReference

gt = ""
trs = arcpy.ListTransformations(in_sr, out_sr)
if trs:
    gt = trs[0]  # ArcGIS-picked best available on your machine

# Use template cell size numbers (meters)
cellx = abs(arcpy.Describe(TEMPLATE).meanCellWidth)
celly = abs(arcpy.Describe(TEMPLATE).meanCellHeight)
cellsize = f"{cellx} {celly}"

# Limit extent to template (important for EU-scale rasters)
ext = arcpy.Describe(TEMPLATE).extent
arcpy.env.extent = ext  # for ProjectRaster output footprint

arcpy.management.ProjectRaster(
    in_raster=MMQ_CLEAN_GCS,
    out_raster=MMQ_PROJ_RAW,
    out_coor_system=out_sr,
    resampling_type="BILINEAR",   # discharge is continuous
    cell_size=cellsize,
    geographic_transform=gt
)

# 3) Snap/align EXACTLY to template grid (origin, rows/cols)
arcpy.env.snapRaster = TEMPLATE
arcpy.env.cellSize   = TEMPLATE
arcpy.env.extent     = TEMPLATE

arcpy.management.Resample(
    in_raster=MMQ_PROJ_RAW,
    out_raster=MMQ_ALIGNED,
    cell_size=cellsize,
    resampling_type="BILINEAR"
)

# 4) Re-clean (optional safety) after resampling
mmq2 = Raster(MMQ_ALIGNED)
mmq2_clean = SetNull((mmq2 < 0), mmq2)
mmq2_clean.save(MMQ_FINAL)

print("DONE:", MMQ_FINAL)
print("CRS:", arcpy.Describe(MMQ_FINAL).spatialReference.name)
print("Extent:", arcpy.Describe(MMQ_FINAL).extent)
print("Cell:", arcpy.Describe(MMQ_FINAL).meanCellWidth, arcpy.Describe(MMQ_FINAL).meanCellHeight)
