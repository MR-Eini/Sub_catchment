#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from pathlib import Path
import arcpy
from arcpy.sa import Raster, ZonalStatistics, SetNull, IsNull, Con

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True

# -----------------------------
# EDIT PATHS
# -----------------------------
BASE = r"E:\pythonProject\Europe\OUT_ARCGIS_90M"
SUBCATCH = os.path.join(BASE, "subcatch_streamlink.tif")

MMQ_PATH = r"E:\pythonProject\Europe\DELTAQ_90M_MMQ_FIX\MMQ_90m_3035_aligned.tif"

OUT_DIR = r"E:\pythonProject\Europe\DELTAQ_90M_MMQ_FIX"
Path(OUT_DIR).mkdir(parents=True, exist_ok=True)

NODATA_VAL = -9999

MMQ_CLEAN0_OUT    = os.path.join(OUT_DIR, "MMQ_clean0.tif")                 # optional intermediate
MMQ_SUBCATCH_OUT  = os.path.join(OUT_DIR, "MMQ_mean_by_subcatch.tif")       # final

# -----------------------------
# ENV: force SUBCATCH grid
# -----------------------------
arcpy.env.snapRaster = SUBCATCH
arcpy.env.cellSize   = SUBCATCH
arcpy.env.extent     = SUBCATCH

# -----------------------------
# 0) Prepare SUBCATCH zones (exclude 0 or NoData zones)
# -----------------------------
sub = Raster(SUBCATCH)
sub_zones = SetNull(IsNull(sub) | (sub <= 0), sub)  # only real catchment IDs remain

# -----------------------------
# 1) FIRST: change MMQ -9999 (and true NoData) to 0
# -----------------------------
mmq = Raster(MMQ_PATH)

mmq_clean0 = Con(IsNull(mmq) | (mmq == NODATA_VAL), 0, mmq)   # all missing -> 0
# safety: if anything still equals -9999 (e.g., carried through), force to 0
mmq_clean0 = Con(mmq_clean0 == NODATA_VAL, 0, mmq_clean0)

# optional save intermediate
if arcpy.Exists(MMQ_CLEAN0_OUT):
    arcpy.management.Delete(MMQ_CLEAN0_OUT)
mmq_clean0.save(MMQ_CLEAN0_OUT)

# -----------------------------
# 2) DO NOT calculate mean where MMQ == 0
#    (exclude zeros from averaging)
# -----------------------------
mmq_for_mean = SetNull(mmq_clean0 == 0, mmq_clean0)  # zeros ignored in MEAN

# -----------------------------
# 3) Mean MMQ per catchment (average of nonzero MMQ pixels inside each zone)
# -----------------------------
mmq_mean = ZonalStatistics(sub_zones, "Value", mmq_for_mean, "MEAN", "DATA")

# -----------------------------
# 4) Fill everywhere with 0 where mean is missing OR outside catchments
# -----------------------------
mmq_mean0 = Con(IsNull(mmq_mean) | (mmq_mean == NODATA_VAL), 0, mmq_mean)
final = Con(IsNull(sub_zones), 0, mmq_mean0)

# final safety pass: no -9999 anywhere
final = Con(final == NODATA_VAL, 0, final)

# -----------------------------
# 5) Save
# -----------------------------
if arcpy.Exists(MMQ_SUBCATCH_OUT):
    arcpy.management.Delete(MMQ_SUBCATCH_OUT)

final.save(MMQ_SUBCATCH_OUT)

print("Saved MMQ_clean0:", MMQ_CLEAN0_OUT)
print("Saved MMQ mean by subcatch (no -9999, no NoData):", MMQ_SUBCATCH_OUT)
