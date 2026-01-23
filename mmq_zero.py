#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from pathlib import Path
import arcpy
from arcpy.sa import Raster, SetNull

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True

# -----------------------------
# EDIT THESE
# -----------------------------
IN_DIR  = r"/Europe/MMQ_clean"
# only process rasters whose filename contains this text (set to "" to process all .tif)
NAME_CONTAINS = "MMQ_clean0"

OUT_DIR = os.path.join(IN_DIR, "MMQ_clean0_0_to_NoData")
Path(OUT_DIR).mkdir(parents=True, exist_ok=True)

# If True, also treat very small values as zero (optional safety)
USE_TOLERANCE = False
ZERO_TOL = 1e-12

# -----------------------------
# PROCESS
# -----------------------------
arcpy.env.workspace = IN_DIR

tifs = []
for root, _, files in os.walk(IN_DIR):
    for fn in files:
        if not fn.lower().endswith(".tif"):
            continue
        if NAME_CONTAINS and (NAME_CONTAINS not in fn):
            continue
        tifs.append(os.path.join(root, fn))

if not tifs:
    raise RuntimeError(f"No matching .tif found in: {IN_DIR}")

print(f"Found {len(tifs)} raster(s). Output -> {OUT_DIR}")

for in_path in tifs:
    base = os.path.splitext(os.path.basename(in_path))[0]
    out_path = os.path.join(OUT_DIR, f"{base}_0NoData.tif")

    if arcpy.Exists(out_path):
        arcpy.management.Delete(out_path)

    r = Raster(in_path)

    if USE_TOLERANCE:
        # Set values with abs(value) <= ZERO_TOL to NoData
        r_out = SetNull(abs(r) <= ZERO_TOL, r)
    else:
        # Set exactly 0 to NoData
        # (SetNull(condition_raster, false_raster, where_clause))
        r_out = SetNull(r, r, "VALUE = 0")

    r_out.save(out_path)

    print("Saved:", out_path)

print("Done.")
