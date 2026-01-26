#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MMQ flux per StreamLink/subcatchment:
  flux = outlet - inlet

Headwater/single-endpoint rule requested:
  If outlet is missing but inlet exists -> set outlet = inlet and inlet = 0.

Also fixes the "CSV has fewer rows" issue by building a master ID table from the
StreamLink raster attribute table, then left-joining inlet/outlet stats onto it.
"""

import os
import arcpy
from arcpy.sa import Raster, SetNull, Con, IsNull, ZonalStatistics, ExtractValuesToPoints

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True

# -----------------------------
# INPUTS
# -----------------------------
MMQ      = r"E:\pythonProject\Europe\MMQ_ALIGNED_25M\MMQ_mamad_3035_25m_clean.tif"
FACC     = r"E:\pythonProject\Europe\OUT_ARCGIS_25M\flowacc_cells.tif"
SLINK    = r"E:\pythonProject\Europe\OUT_ARCGIS_25M\streamlink.tif"
SUBCATCH = r"E:\pythonProject\Europe\OUT_ARCGIS_25M\subcatch_streamlink.tif"  # not used here, but kept for completeness

OUT_DIR = r"E:\pythonProject\Europe\MMQ_FLUX_25M"
os.makedirs(OUT_DIR, exist_ok=True)

# If MMQ_out == 0 should be considered "missing outlet", set True
OUTLET_ZERO_COUNTS_AS_MISSING = True

# -----------------------------
# WORKSPACE
# -----------------------------
GDB = os.path.join(OUT_DIR, "mmq_flux.gdb")
if not arcpy.Exists(GDB):
    arcpy.management.CreateFileGDB(OUT_DIR, "mmq_flux.gdb")

arcpy.env.workspace = GDB
arcpy.env.scratchWorkspace = GDB

# Align to SLINK grid
arcpy.env.snapRaster = SLINK
arcpy.env.cellSize   = SLINK
arcpy.env.extent     = SLINK

# -----------------------------
# OUTPUTS
# -----------------------------
mmq_clean_ras   = os.path.join(GDB, "mmq_clean")

zmax_facc_ras   = os.path.join(GDB, "zmax_facc")
zmin_facc_ras   = os.path.join(GDB, "zmin_facc")

outlet_id_ras   = os.path.join(GDB, "outlet_id")
inlet_id_ras    = os.path.join(GDB, "inlet_id")

outlet_pts      = os.path.join(GDB, "outlet_pts")
inlet_pts       = os.path.join(GDB, "inlet_pts")

outlet_pts_ev   = os.path.join(GDB, "outlet_pts_mmq")
inlet_pts_ev    = os.path.join(GDB, "inlet_pts_mmq")

outlet_stats    = os.path.join(GDB, "outlet_stats")
inlet_stats     = os.path.join(GDB, "inlet_stats")

all_ids_tbl     = os.path.join(GDB, "all_slink_ids")
flux_table      = os.path.join(GDB, "mmq_flux_by_slink")
flux_csv        = os.path.join(OUT_DIR, "mmq_flux_by_slink.csv")


def ensure_aligned(name, ras_path):
    d = arcpy.Describe(ras_path)
    print(f"{name}: {ras_path}")
    print(f"  CRS: {d.spatialReference.name}")
    print(f"  Cell: {d.meanCellWidth} {d.meanCellHeight}")
    print(f"  Ext:  {d.extent.XMin:.3f} {d.extent.YMin:.3f} {d.extent.XMax:.3f} {d.extent.YMax:.3f}")


def rename_gridcode_to_slinkid(fc):
    fields = arcpy.ListFields(fc)
    gc = [f.name for f in fields if f.name.lower() == "grid_code"]
    if gc:
        arcpy.management.AlterField(fc, gc[0], "SLINK_ID", "SLINK_ID")


def main():
    ensure_aligned("MMQ", MMQ)
    ensure_aligned("FACC", FACC)
    ensure_aligned("SLINK", SLINK)

    # 0) Clean MMQ
    print("\n0) Cleaning MMQ (-9999 or <0 -> NoData)...")
    mmq_r = Raster(MMQ)
    mmq_clean = SetNull((mmq_r == -9999) | (mmq_r < 0), mmq_r)
    mmq_clean.save(mmq_clean_ras)

    # 1) Zonal MAX/MIN of FACC per StreamLink ID
    print("\n1) ZonalStatistics on FACC by StreamLink (MAXIMUM/MINIMUM)...")
    zmax = ZonalStatistics(SLINK, "Value", Raster(FACC), "MAXIMUM", "DATA")
    zmin = ZonalStatistics(SLINK, "Value", Raster(FACC), "MINIMUM", "DATA")
    zmax.save(zmax_facc_ras)
    zmin.save(zmin_facc_ras)

    # 2) Endpoint rasters (IDs at endpoint cells)
    print("\n2) Building inlet/outlet ID rasters...")
    sl = Raster(SLINK)
    fa = Raster(FACC)

    is_stream = ~IsNull(sl)
    outlet_mask = (fa == Raster(zmax_facc_ras)) & is_stream
    inlet_mask  = (fa == Raster(zmin_facc_ras)) & is_stream

    outlet_id = Con(outlet_mask, sl)
    inlet_id  = Con(inlet_mask,  sl)

    outlet_id.save(outlet_id_ras)
    inlet_id.save(inlet_id_ras)

    # 3) RasterToPoint endpoints
    print("\n3) RasterToPoint endpoints...")
    if arcpy.Exists(outlet_pts): arcpy.management.Delete(outlet_pts)
    if arcpy.Exists(inlet_pts):  arcpy.management.Delete(inlet_pts)

    arcpy.conversion.RasterToPoint(outlet_id_ras, outlet_pts, "VALUE")
    arcpy.conversion.RasterToPoint(inlet_id_ras,  inlet_pts,  "VALUE")

    rename_gridcode_to_slinkid(outlet_pts)
    rename_gridcode_to_slinkid(inlet_pts)

    # 4) Sample MMQ at endpoints
    print("\n4) ExtractValuesToPoints (MMQ at inlet/outlet points)...")
    if arcpy.Exists(outlet_pts_ev): arcpy.management.Delete(outlet_pts_ev)
    if arcpy.Exists(inlet_pts_ev):  arcpy.management.Delete(inlet_pts_ev)

    ExtractValuesToPoints(outlet_pts, mmq_clean_ras, outlet_pts_ev, "NONE", "VALUE_ONLY")
    ExtractValuesToPoints(inlet_pts,  mmq_clean_ras, inlet_pts_ev,  "NONE", "VALUE_ONLY")

    mmq_field = "RASTERVALU"

    # 5) Stats per SLINK_ID (mean handles ties)
    print("\n5) Statistics by SLINK_ID (mean MMQ at endpoints)...")
    if arcpy.Exists(outlet_stats): arcpy.management.Delete(outlet_stats)
    if arcpy.Exists(inlet_stats):  arcpy.management.Delete(inlet_stats)

    arcpy.analysis.Statistics(outlet_pts_ev, outlet_stats, [[mmq_field, "MEAN"]], "SLINK_ID")
    arcpy.analysis.Statistics(inlet_pts_ev,  inlet_stats,  [[mmq_field, "MEAN"]], "SLINK_ID")

    out_mean = [f.name for f in arcpy.ListFields(outlet_stats) if f.name.upper().startswith("MEAN_")][0]
    in_mean  = [f.name for f in arcpy.ListFields(inlet_stats)  if f.name.upper().startswith("MEAN_")][0]

    arcpy.management.AlterField(outlet_stats, out_mean, "MMQ_out", "MMQ_out")
    arcpy.management.AlterField(inlet_stats,  in_mean,  "MMQ_in",  "MMQ_in")

    # 6) Build master ID table from StreamLink RAT so CSV has ALL IDs
    print("\n6) Building master ID table from StreamLink raster attribute table...")
    arcpy.management.BuildRasterAttributeTable(SLINK, "Overwrite")

    if arcpy.Exists(all_ids_tbl): arcpy.management.Delete(all_ids_tbl)
    arcpy.management.CopyRows(SLINK, all_ids_tbl)  # copies RAT (Value, Count, ...)

    # Rename Value -> SLINK_ID
    if "Value" not in [f.name for f in arcpy.ListFields(all_ids_tbl)]:
        raise RuntimeError("Could not find 'Value' field in StreamLink raster attribute table.")
    arcpy.management.AlterField(all_ids_tbl, "Value", "SLINK_ID", "SLINK_ID")

    # 7) Join outlet/inlet onto master table, then apply requested rule
    print("\n7) Joining inlet/outlet and computing MMQ_flux with headwater rule...")
    if arcpy.Exists(flux_table): arcpy.management.Delete(flux_table)
    arcpy.management.CopyRows(all_ids_tbl, flux_table)

    arcpy.management.JoinField(flux_table, "SLINK_ID", outlet_stats, "SLINK_ID", ["MMQ_out"])
    arcpy.management.JoinField(flux_table, "SLINK_ID", inlet_stats,  "SLINK_ID", ["MMQ_in"])

    # Add corrected fields
    for fld, ftype, flen in [
        ("MMQ_out_fix", "DOUBLE", None),
        ("MMQ_in_fix",  "DOUBLE", None),
        ("MMQ_flux",    "DOUBLE", None),
        ("FLAG",        "TEXT",   50),
    ]:
        if fld not in [f.name for f in arcpy.ListFields(flux_table)]:
            if ftype == "TEXT":
                arcpy.management.AddField(flux_table, fld, ftype, field_length=flen)
            else:
                arcpy.management.AddField(flux_table, fld, ftype)

    # Apply rule:
    # - if no outlet but inlet exists -> outlet=inlet ; inlet=0
    # - if no inlet but outlet exists -> inlet=0
    # - otherwise -> keep both
    with arcpy.da.UpdateCursor(
        flux_table,
        ["MMQ_out", "MMQ_in", "MMQ_out_fix", "MMQ_in_fix", "MMQ_flux", "FLAG"]
    ) as cur:
        for outv, inv, out_fix, in_fix, flx, flag in cur:

            outlet_missing = (outv is None) or (OUTLET_ZERO_COUNTS_AS_MISSING and outv == 0)
            inlet_missing  = (inv is None)

            if outlet_missing and (not inlet_missing):
                # requested: use inlet as outlet
                out_use = float(inv)
                in_use  = 0.0
                flag_use = "OUT_MISSING_USE_IN_AS_OUT"
            else:
                out_use = None if outv is None else float(outv)
                in_use  = 0.0 if inv is None else float(inv)
                flag_use = "IN_MISSING_ASSUMED0" if (not outlet_missing and inlet_missing) else ""

            flux_use = None if out_use is None else (out_use - in_use)

            cur.updateRow((outv, inv, out_use, in_use, flux_use, flag_use))

    # Export CSV
    if os.path.exists(flux_csv):
        try:
            os.remove(flux_csv)
        except Exception:
            pass
    arcpy.conversion.TableToTable(flux_table, OUT_DIR, os.path.basename(flux_csv))
    print(f"\nDONE:\n  Table: {flux_table}\n  CSV  : {flux_csv}")


if __name__ == "__main__":
    main()
