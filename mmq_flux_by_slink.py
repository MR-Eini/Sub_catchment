#!/usr/bin/env python3
import os
import arcpy
from arcpy.sa import (
    Raster, SetNull, Con, IsNull, ZonalStatistics,
    ExtractValuesToPoints
)

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True

# -----------------------------
# INPUTS
# -----------------------------
MMQ     = r"E:\pythonProject\Europe\MMQ_ALIGNED_25M\MMQ_mamad_3035_25m_clean.tif"
FDIR    = r"E:\pythonProject\Europe\OUT_ARCGIS_25M\flowdir_d8.tif"          # not used directly, but good for snapping
FACC    = r"E:\pythonProject\Europe\OUT_ARCGIS_25M\flowacc_cells.tif"
SLINK   = r"E:\pythonProject\Europe\OUT_ARCGIS_25M\streamlink.tif"
SUBCATCH= r"E:\pythonProject\Europe\OUT_ARCGIS_25M\subcatch_streamlink.tif" # same IDs as SLINK

OUT_DIR = r"E:\pythonProject\Europe\MMQ_FLUX_25M"
os.makedirs(OUT_DIR, exist_ok=True)

# Use a FileGDB for speed / field name safety
GDB = os.path.join(OUT_DIR, "mmq_flux.gdb")
if not arcpy.Exists(GDB):
    arcpy.management.CreateFileGDB(OUT_DIR, "mmq_flux.gdb")

# -----------------------------
# ENVIRONMENT (important!)
# -----------------------------
arcpy.env.snapRaster = SLINK
arcpy.env.cellSize   = SLINK
arcpy.env.extent     = SLINK
arcpy.env.workspace  = GDB
arcpy.env.scratchWorkspace = GDB

# Outputs
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

flux_table      = os.path.join(GDB, "mmq_flux_by_slink")
flux_csv        = os.path.join(OUT_DIR, "mmq_flux_by_slink.csv")

# Optional flux raster (constant per subcatch)
flux_reclass_tbl= os.path.join(GDB, "reclass_tbl")
flux_raster     = os.path.join(OUT_DIR, "MMQ_flux_subcatch.tif")

def ensure_aligned(name, ras_path):
    """Basic alignment sanity (CRS, cellsize); if mismatch, project/resample outside this script."""
    d = arcpy.Describe(ras_path)
    # You can add stricter checks here if needed.
    print(f"{name}: {ras_path}")
    print(f"  CRS: {d.spatialReference.name}")
    print(f"  Ext: {d.extent.XMin:.3f} {d.extent.YMin:.3f} {d.extent.XMax:.3f} {d.extent.YMax:.3f}")

def main():
    ensure_aligned("MMQ", MMQ)
    ensure_aligned("FACC", FACC)
    ensure_aligned("SLINK", SLINK)

    # 0) Clean MMQ: set -9999 and negatives to NoData (adjust if your NoData differs)
    print("\n0) Cleaning MMQ (-9999 or <0 -> NoData)...")
    mmq_r = Raster(MMQ)
    mmq_clean = SetNull((mmq_r == -9999) | (mmq_r < 0), mmq_r)
    mmq_clean.save(mmq_clean_ras)

    # 1) Zonal MAX/MIN of flowacc along each StreamLink ID
    #    (SLINK raster typically exists only on stream cells; this is what we want)
    print("\n1) ZonalStatistics on FACC by StreamLink (MAX/MIN)...")
    zmax = ZonalStatistics(SLINK, "Value", Raster(FACC), "MAXIMUM", "DATA")
    zmin = ZonalStatistics(SLINK, "Value", Raster(FACC), "MINIMUM", "DATA")
    zmax.save(zmax_facc_ras)
    zmin.save(zmin_facc_ras)

    # 2) Endpoint masks: cells where FACC equals zonal max/min inside each link
    print("\n2) Building inlet/outlet ID rasters (one cell per endpoint, per link; ties possible)...")
    sl = Raster(SLINK)
    fa = Raster(FACC)

    is_stream = ~IsNull(sl)

    outlet_mask = (fa == Raster(zmax_facc_ras)) & is_stream
    inlet_mask  = (fa == Raster(zmin_facc_ras)) & is_stream

    outlet_id = Con(outlet_mask, sl)  # cell value = StreamLink ID at outlet cell(s)
    inlet_id  = Con(inlet_mask,  sl)  # cell value = StreamLink ID at inlet cell(s)

    outlet_id.save(outlet_id_ras)
    inlet_id.save(inlet_id_ras)

    # 3) Raster endpoints -> points (grid_code = StreamLink ID)
    print("\n3) RasterToPoint endpoints...")
    arcpy.conversion.RasterToPoint(outlet_id_ras, outlet_pts, "VALUE")
    arcpy.conversion.RasterToPoint(inlet_id_ras,  inlet_pts,  "VALUE")

    # Rename grid_code to SLINK_ID (field names in GDB are fine)
    for fc in [outlet_pts, inlet_pts]:
        if "grid_code" in [f.name.lower() for f in arcpy.ListFields(fc)]:
            # ArcGIS keeps exact case; find actual field name
            gc = [f.name for f in arcpy.ListFields(fc) if f.name.lower() == "grid_code"][0]
            arcpy.management.AlterField(fc, gc, "SLINK_ID", "SLINK_ID")

    # 4) Sample MMQ at those endpoints
    print("\n4) ExtractValuesToPoints (MMQ at inlet/outlet points)...")
    ExtractValuesToPoints(outlet_pts, mmq_clean_ras, outlet_pts_ev, "NONE", "VALUE_ONLY")
    ExtractValuesToPoints(inlet_pts,  mmq_clean_ras, inlet_pts_ev,  "NONE", "VALUE_ONLY")

    # ExtractValuesToPoints adds field "RASTERVALU" (MMQ)
    mmq_field = "RASTERVALU"

    # 5) If ties created multiple points per SLINK_ID, collapse by mean (safe)
    print("\n5) Statistics by SLINK_ID (mean MMQ at endpoints, handling ties)...")
    arcpy.analysis.Statistics(outlet_pts_ev, outlet_stats, [[mmq_field, "MEAN"]], "SLINK_ID")
    arcpy.analysis.Statistics(inlet_pts_ev,  inlet_stats,  [[mmq_field, "MEAN"]], "SLINK_ID")

    # Rename MEAN_RASTERVALU to MMQ_out / MMQ_in
    out_mean = [f.name for f in arcpy.ListFields(outlet_stats) if f.name.upper().startswith("MEAN_")][0]
    in_mean  = [f.name for f in arcpy.ListFields(inlet_stats)  if f.name.upper().startswith("MEAN_")][0]

    arcpy.management.AlterField(outlet_stats, out_mean, "MMQ_out", "MMQ_out")
    arcpy.management.AlterField(inlet_stats,  in_mean,  "MMQ_in",  "MMQ_in")

    # 6) Join inlet to outlet and compute flux
    print("\n6) Joining inlet/outlet and computing MMQ_flux = out - in ...")
    arcpy.management.CopyRows(outlet_stats, flux_table)
    arcpy.management.JoinField(flux_table, "SLINK_ID", inlet_stats, "SLINK_ID", ["MMQ_in"])

    # Add flux field
    if "MMQ_flux" not in [f.name for f in arcpy.ListFields(flux_table)]:
        arcpy.management.AddField(flux_table, "MMQ_flux", "DOUBLE")

    # Treat missing inlet as 0 (typical for headwaters OR NoData inlet sampling)
    arcpy.management.CalculateField(
        flux_table,
        "MMQ_flux",
        "!MMQ_out! - (0 if !MMQ_in! is None else !MMQ_in!)",
        "PYTHON3"
    )

    # Export CSV
    arcpy.conversion.TableToTable(flux_table, OUT_DIR, os.path.basename(flux_csv))
    print(f"\nTable written:\n  {flux_table}\n  {flux_csv}")

    # 7) OPTIONAL: Creating flux raster per subcatch (ReclassByTable)...
    print("\n7) OPTIONAL: Creating flux raster per subcatch (ReclassByTable)...")

    # Create a 3-field table: FROM_VAL, TO_VAL, MMQ_flux  (FROM=TO=ID)
    safe_tbl_name = "reclass_tbl"
    flux_reclass_tbl = os.path.join(GDB, safe_tbl_name)

    if arcpy.Exists(flux_reclass_tbl):
        arcpy.management.Delete(flux_reclass_tbl)

    arcpy.management.CreateTable(GDB, safe_tbl_name)

    arcpy.management.AddField(flux_reclass_tbl, "FROM_VAL", "LONG")
    arcpy.management.AddField(flux_reclass_tbl, "TO_VAL", "LONG")
    arcpy.management.AddField(flux_reclass_tbl, "MMQ_flux", "DOUBLE")

    with arcpy.da.InsertCursor(flux_reclass_tbl, ["FROM_VAL", "TO_VAL", "MMQ_flux"]) as icur:
        with arcpy.da.SearchCursor(flux_table, ["SLINK_ID", "MMQ_flux"]) as scur:
            for sid, flx in scur:
                sid_i = int(sid)
                icur.insertRow((sid_i, sid_i, float(flx) if flx is not None else None))

    flux_sub = arcpy.sa.ReclassByTable(
        in_raster=SUBCATCH,
        in_remap_table=flux_reclass_tbl,
        from_value_field="FROM_VAL",
        to_value_field="TO_VAL",
        output_value_field="MMQ_flux",
        missing_values="NODATA"
    )
    flux_sub.save(flux_raster)
    print(f"Flux raster written:\n  {flux_raster}")

if __name__ == "__main__":
    main()
