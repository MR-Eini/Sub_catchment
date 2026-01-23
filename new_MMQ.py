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

MMQ_CLEAN0_OUT   = os.path.join(OUT_DIR, "MMQ_clean0.tif")            # intermediate (NoData/neg -> 0)
MMQ_SUBCATCH_OUT = os.path.join(OUT_DIR, "MMQ_mean_by_subcatch.tif")  # final (mean per zone, filled with 0)

# -----------------------------
# ENV: force SUBCATCH grid
# -----------------------------
arcpy.env.snapRaster = SUBCATCH
arcpy.env.cellSize   = SUBCATCH
arcpy.env.extent     = SUBCATCH
arcpy.env.mask       = SUBCATCH  # optional but helpful

def rp(path_or_ras, prop):
    return arcpy.management.GetRasterProperties(path_or_ras, prop).getOutput(0)

# -----------------------------
# 0) Prepare SUBCATCH zones (exclude 0/NoData)
# -----------------------------
sub = Raster(SUBCATCH)
sub_zones = SetNull(IsNull(sub) | (sub <= 0), sub)

# Diagnostics: CRS + extent
d_sub = arcpy.Describe(SUBCATCH)
d_mmq = arcpy.Describe(MMQ_PATH)
print("SUBCATCH CRS:", d_sub.spatialReference.name)
print("MMQ CRS     :", d_mmq.spatialReference.name)
print("SUBCATCH extent:", d_sub.extent)
print("MMQ extent     :", d_mmq.extent)

# -----------------------------
# 1) Save MMQ with NoData/negative -> 0 (robust)
#    NOTE: your resampled NoData is ~ -9998.999, so we treat ANY negative as missing.
# -----------------------------
mmq = Raster(MMQ_PATH)

mmq_clean0 = Con(IsNull(mmq) | (mmq < 0), 0, mmq)
mmq_clean0 = Con(mmq_clean0 < 0, 0, mmq_clean0)  # safety

if arcpy.Exists(MMQ_CLEAN0_OUT):
    arcpy.management.Delete(MMQ_CLEAN0_OUT)
mmq_clean0.save(MMQ_CLEAN0_OUT)

print("MMQ_clean0 MIN/MAX:", rp(MMQ_CLEAN0_OUT, "MINIMUM"), rp(MMQ_CLEAN0_OUT, "MAXIMUM"))

# -----------------------------
# 2) For MEAN: ignore all <= 0 (and any NoData)
#    Compute from ORIGINAL mmq to avoid any pixel-type rounding in saved intermediate.
# -----------------------------
mmq_for_mean = SetNull(IsNull(mmq) | (mmq <= 0), mmq)

# Optional diagnostics (may fail if no stats; fine)
try:
    print("mmq_for_mean MIN/MAX:", rp(mmq_for_mean, "MINIMUM"), rp(mmq_for_mean, "MAXIMUM"))
except Exception as e:
    print("mmq_for_mean MIN/MAX: could not read raster properties:", e)

# -----------------------------
# 3) Mean MMQ per catchment (only from positive MMQ pixels inside each zone)
# -----------------------------
mmq_mean = ZonalStatistics(sub_zones, "Value", mmq_for_mean, "MEAN", "DATA")

# -----------------------------
# 4) Fill with 0 where mean is missing OR outside catchments
# -----------------------------
final = Con(IsNull(sub_zones) | IsNull(mmq_mean), 0, mmq_mean)

# Final safety: no negatives anywhere
final = Con(final < 0, 0, final)

# -----------------------------
# 5) Save
# -----------------------------
if arcpy.Exists(MMQ_SUBCATCH_OUT):
    arcpy.management.Delete(MMQ_SUBCATCH_OUT)

final.save(MMQ_SUBCATCH_OUT)

print("Saved MMQ_clean0:", MMQ_CLEAN0_OUT)
print("Saved MMQ mean by subcatch (no NoData, no negatives):", MMQ_SUBCATCH_OUT)
print("FINAL MIN/MAX:", rp(MMQ_SUBCATCH_OUT, "MINIMUM"), rp(MMQ_SUBCATCH_OUT, "MAXIMUM"))
