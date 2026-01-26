#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute Qmin flux per StreamLink/subcatchment:
  Qmin_flux = Qmin_out - Qmin_in

Headwater rule (as requested):
  If no outlet but inlet exists -> outlet = inlet ; inlet = 0.

This script:
1) (If needed) projects + resamples Qmin to match streamlink grid (EPSG:3035, 90m, snap/origin)
2) extracts Qmin at inlet/outlet (by FACC min/max within each StreamLink)
3) builds a master ID table from StreamLink RAT (so ALL IDs are present)
4) computes Qmin_flux with the headwater rule
5) exports CSV using Python (no lock problems)
6) optional: creates dissolved polygons + joined shapefile + Qmin_flux raster
"""

import os
import csv
import arcpy
from arcpy.sa import Raster, SetNull, Con, IsNull, ZonalStatistics, ExtractValuesToPoints

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True

# -----------------------------
# INPUTS
# -----------------------------
QMIN_IN   = r"E:\pythonProject\Europe\Qmin_mamad.tif"
FACC      = r"E:\pythonProject\Europe\OUT_ARCGIS_25M\flowacc_cells.tif"
SLINK     = r"E:\pythonProject\Europe\OUT_ARCGIS_25M\streamlink.tif"
SUBCATCH  = r"E:\pythonProject\Europe\OUT_ARCGIS_25M\subcatch_streamlink.tif"

OUT_DIR   = r"E:\pythonProject\Europe\QMIN_FLUX_25M"
os.makedirs(OUT_DIR, exist_ok=True)

# -----------------------------
# OPTIONS
# -----------------------------
TREAT_ZERO_AS_MISSING = False          # Qmin=0 can be valid; keep False unless 0 means NoData in your raster
OUTLET_ZERO_COUNTS_AS_MISSING = False  # set True only if you want to treat outlet=0 as "missing outlet"
MAKE_SHP   = True
MAKE_RASTER= True

# -----------------------------
# WORKSPACE
# -----------------------------
GDB = os.path.join(OUT_DIR, "qmin_flux.gdb")
if not arcpy.Exists(GDB):
    arcpy.management.CreateFileGDB(OUT_DIR, "qmin_flux.gdb")
arcpy.env.workspace = GDB
arcpy.env.scratchWorkspace = GDB

# template alignment
arcpy.env.snapRaster = SLINK
arcpy.env.cellSize   = SLINK
arcpy.env.extent     = SLINK

# -----------------------------
# OUTPUTS
# -----------------------------
ALN_DIR        = os.path.join(OUT_DIR, "aligned")
os.makedirs(ALN_DIR, exist_ok=True)

QMIN_CLEAN_SRC = os.path.join(ALN_DIR, "Qmin_clean_src.tif")
QMIN_PROJ_RAW  = os.path.join(ALN_DIR, "Qmin_3035_raw.tif")
QMIN_ALIGNED   = os.path.join(ALN_DIR, "Qmin_3035_aligned.tif")

qmin_clean_ras = os.path.join(GDB, "qmin_clean")

zmax_facc_ras  = os.path.join(GDB, "zmax_facc")
zmin_facc_ras  = os.path.join(GDB, "zmin_facc")

outlet_id_ras  = os.path.join(GDB, "outlet_id")
inlet_id_ras   = os.path.join(GDB, "inlet_id")

outlet_pts     = os.path.join(GDB, "outlet_pts")
inlet_pts      = os.path.join(GDB, "inlet_pts")

outlet_pts_ev  = os.path.join(GDB, "outlet_pts_qmin")
inlet_pts_ev   = os.path.join(GDB, "inlet_pts_qmin")

outlet_stats   = os.path.join(GDB, "outlet_stats")
inlet_stats    = os.path.join(GDB, "inlet_stats")

all_ids_tbl    = os.path.join(GDB, "all_slink_ids")
flux_table     = os.path.join(GDB, "qmin_flux_by_slink")

CSV_OUT        = os.path.join(OUT_DIR, "qmin_flux_by_slink.csv")

poly_raw       = os.path.join(GDB, "subcatch_poly_raw")
poly_diss      = os.path.join(GDB, "subcatch_poly_diss")
poly_join      = os.path.join(GDB, "subcatch_poly_join")

SHP_OUT        = os.path.join(OUT_DIR, "subcatch_with_qmin_flux.shp")
FLUX_RAS_OUT   = os.path.join(OUT_DIR, "Qmin_flux_subcatch.tif")


def desc_ras(path):
    d = arcpy.Describe(path)
    return d.spatialReference.name, float(d.meanCellWidth), float(d.meanCellHeight), d.extent


def ensure_aligned(name, ras_path):
    crs, cx, cy, e = desc_ras(ras_path)
    print(f"{name}: {ras_path}")
    print(f"  CRS: {crs}")
    print(f"  Cell: {cx} {cy}")
    print(f"  Ext:  {e.XMin:.3f} {e.YMin:.3f} {e.XMax:.3f} {e.YMax:.3f}")


def is_same_grid(a, b):
    crs_a, cx_a, cy_a, e_a = desc_ras(a)
    crs_b, cx_b, cy_b, e_b = desc_ras(b)
    if crs_a != crs_b:
        return False
    if abs(cx_a - cx_b) > 1e-6 or abs(cy_a - cy_b) > 1e-6:
        return False
    # extent equality is tricky; snapping/extent env will handle footprint
    return True


def align_to_template(in_ras, template_ras):
    """
    Project + resample (snap) to template grid if needed.
    Returns aligned raster path.
    """
    if is_same_grid(in_ras, template_ras):
        return in_ras

    in_sr  = arcpy.Describe(in_ras).spatialReference
    out_sr = arcpy.Describe(template_ras).spatialReference

    # 1) Clean in source CRS (avoid -9999 interpolation)
    r = Raster(in_ras)
    if TREAT_ZERO_AS_MISSING:
        r_clean = SetNull((r == -9999) | (r < 0) | (r == 0), r)
    else:
        r_clean = SetNull((r == -9999) | (r < 0), r)
    r_clean.save(QMIN_CLEAN_SRC)

    # choose a geographic transformation if available
    gt = ""
    trs = arcpy.ListTransformations(in_sr, out_sr)
    if trs:
        gt = trs[0]

    # template cell size
    dtemp = arcpy.Describe(template_ras)
    cellx = abs(float(dtemp.meanCellWidth))
    celly = abs(float(dtemp.meanCellHeight))
    cellsize = f"{cellx} {celly}"

    # project (continuous -> bilinear)
    arcpy.env.extent = dtemp.extent
    arcpy.management.ProjectRaster(
        in_raster=QMIN_CLEAN_SRC,
        out_raster=QMIN_PROJ_RAW,
        out_coor_system=out_sr,
        resampling_type="BILINEAR",
        cell_size=cellsize,
        geographic_transform=gt
    )

    # resample with snap to ensure identical origin/grid
    arcpy.env.snapRaster = template_ras
    arcpy.env.cellSize   = template_ras
    arcpy.env.extent     = template_ras

    arcpy.management.Resample(
        in_raster=QMIN_PROJ_RAW,
        out_raster=QMIN_ALIGNED,
        cell_size=cellsize,
        resampling_type="BILINEAR"
    )

    return QMIN_ALIGNED


def rename_gridcode_to_slinkid(fc):
    fields = arcpy.ListFields(fc)
    match = [f.name for f in fields if f.name.lower() == "grid_code"]
    if match:
        arcpy.management.AlterField(fc, match[0], "SLINK_ID", "SLINK_ID")


def safe_csv_path(path):
    if not os.path.exists(path):
        return path
    root, ext = os.path.splitext(path)
    for i in range(1, 1000):
        cand = f"{root}_{i}{ext}"
        if not os.path.exists(cand):
            return cand
    raise RuntimeError("No free CSV filename found.")


def export_table_to_csv(table, out_csv, preferred=None):
    all_fields = [f.name for f in arcpy.ListFields(table) if f.type not in ("Geometry", "Blob", "Raster")]
    if preferred:
        use = [f for f in preferred if f in all_fields]
        if not use:
            use = all_fields
    else:
        use = all_fields

    out_csv = safe_csv_path(out_csv)
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(use)
        with arcpy.da.SearchCursor(table, use) as cur:
            for row in cur:
                w.writerow(row)
    print("CSV:", out_csv)


def main():
    ensure_aligned("QMIN_IN", QMIN_IN)
    ensure_aligned("FACC", FACC)
    ensure_aligned("SLINK", SLINK)

    # A) Align Qmin to template grid if needed
    print("\nA) Aligning Qmin to StreamLink grid (if needed)...")
    QMIN = align_to_template(QMIN_IN, SLINK)
    ensure_aligned("QMIN_USED", QMIN)

    # 0) Clean Qmin again in aligned grid (cheap + consistent)
    print("\n0) Cleaning Qmin in aligned grid...")
    q = Raster(QMIN)
    if TREAT_ZERO_AS_MISSING:
        q_clean = SetNull((q == -9999) | (q < 0) | (q == 0), q)
    else:
        q_clean = SetNull((q == -9999) | (q < 0), q)
    q_clean.save(qmin_clean_ras)

    # 1) Zonal MAX/MIN of FACC per StreamLink
    print("\n1) ZonalStatistics on FACC by StreamLink (MAXIMUM/MINIMUM)...")
    zmax = ZonalStatistics(SLINK, "Value", Raster(FACC), "MAXIMUM", "DATA")
    zmin = ZonalStatistics(SLINK, "Value", Raster(FACC), "MINIMUM", "DATA")
    zmax.save(zmax_facc_ras)
    zmin.save(zmin_facc_ras)

    # 2) Endpoint rasters
    print("\n2) Building inlet/outlet ID rasters...")
    sl = Raster(SLINK)
    fa = Raster(FACC)
    is_stream = ~IsNull(sl)

    outlet_mask = (fa == Raster(zmax_facc_ras)) & is_stream
    inlet_mask  = (fa == Raster(zmin_facc_ras)) & is_stream

    Con(outlet_mask, sl).save(outlet_id_ras)
    Con(inlet_mask,  sl).save(inlet_id_ras)

    # 3) RasterToPoint endpoints
    print("\n3) RasterToPoint endpoints...")
    if arcpy.Exists(outlet_pts): arcpy.management.Delete(outlet_pts)
    if arcpy.Exists(inlet_pts):  arcpy.management.Delete(inlet_pts)
    arcpy.conversion.RasterToPoint(outlet_id_ras, outlet_pts, "VALUE")
    arcpy.conversion.RasterToPoint(inlet_id_ras,  inlet_pts,  "VALUE")
    rename_gridcode_to_slinkid(outlet_pts)
    rename_gridcode_to_slinkid(inlet_pts)

    # 4) Sample Qmin at endpoints
    print("\n4) ExtractValuesToPoints (Qmin at endpoints)...")
    if arcpy.Exists(outlet_pts_ev): arcpy.management.Delete(outlet_pts_ev)
    if arcpy.Exists(inlet_pts_ev):  arcpy.management.Delete(inlet_pts_ev)
    ExtractValuesToPoints(outlet_pts, qmin_clean_ras, outlet_pts_ev, "NONE", "VALUE_ONLY")
    ExtractValuesToPoints(inlet_pts,  qmin_clean_ras, inlet_pts_ev,  "NONE", "VALUE_ONLY")
    val_field = "RASTERVALU"

    # 5) Mean per SLINK_ID (ties)
    print("\n5) Statistics by SLINK_ID (mean at endpoints)...")
    if arcpy.Exists(outlet_stats): arcpy.management.Delete(outlet_stats)
    if arcpy.Exists(inlet_stats):  arcpy.management.Delete(inlet_stats)
    arcpy.analysis.Statistics(outlet_pts_ev, outlet_stats, [[val_field, "MEAN"]], "SLINK_ID")
    arcpy.analysis.Statistics(inlet_pts_ev,  inlet_stats,  [[val_field, "MEAN"]], "SLINK_ID")

    out_mean = [f.name for f in arcpy.ListFields(outlet_stats) if f.name.upper().startswith("MEAN_")][0]
    in_mean  = [f.name for f in arcpy.ListFields(inlet_stats)  if f.name.upper().startswith("MEAN_")][0]
    arcpy.management.AlterField(outlet_stats, out_mean, "Qmin_out", "Qmin_out")
    arcpy.management.AlterField(inlet_stats,  in_mean,  "Qmin_in",  "Qmin_in")

    # 6) Master ID table (ALL IDs)
    print("\n6) Building master ID table from StreamLink RAT...")
    arcpy.management.BuildRasterAttributeTable(SLINK, "Overwrite")
    if arcpy.Exists(all_ids_tbl): arcpy.management.Delete(all_ids_tbl)
    arcpy.management.CopyRows(SLINK, all_ids_tbl)
    if "Value" not in [f.name for f in arcpy.ListFields(all_ids_tbl)]:
        raise RuntimeError("Could not find 'Value' in StreamLink RAT.")
    arcpy.management.AlterField(all_ids_tbl, "Value", "SLINK_ID", "SLINK_ID")

    # 7) Join + headwater rule + flux
    print("\n7) Joining and computing Qmin_flux with rule: if no outlet but inlet -> outlet=inlet, inlet=0")
    if arcpy.Exists(flux_table): arcpy.management.Delete(flux_table)
    arcpy.management.CopyRows(all_ids_tbl, flux_table)
    arcpy.management.JoinField(flux_table, "SLINK_ID", outlet_stats, "SLINK_ID", ["Qmin_out"])
    arcpy.management.JoinField(flux_table, "SLINK_ID", inlet_stats,  "SLINK_ID", ["Qmin_in"])

    for fld, ftype, flen in [
        ("Qmin_out_fix", "DOUBLE", None),
        ("Qmin_in_fix",  "DOUBLE", None),
        ("Qmin_flux",    "DOUBLE", None),
        ("FLAG",         "TEXT",   50),
    ]:
        if fld not in [f.name for f in arcpy.ListFields(flux_table)]:
            if ftype == "TEXT":
                arcpy.management.AddField(flux_table, fld, ftype, field_length=flen)
            else:
                arcpy.management.AddField(flux_table, fld, ftype)

    with arcpy.da.UpdateCursor(
        flux_table,
        ["Qmin_out", "Qmin_in", "Qmin_out_fix", "Qmin_in_fix", "Qmin_flux", "FLAG"]
    ) as cur:
        for outv, inv, out_fix, in_fix, flx, flag in cur:

            outlet_missing = (outv is None) or (OUTLET_ZERO_COUNTS_AS_MISSING and outv == 0)
            inlet_missing  = (inv is None)

            if outlet_missing and (not inlet_missing):
                out_use = float(inv)   # requested
                in_use  = 0.0
                flag_use = "OUT_MISSING_USE_IN_AS_OUT"
            else:
                out_use = None if outv is None else float(outv)
                in_use  = 0.0 if inv is None else float(inv)
                flag_use = "IN_MISSING_ASSUMED0" if (not outlet_missing and inlet_missing) else ""

            flux_use = None if out_use is None else (out_use - in_use)
            cur.updateRow((outv, inv, out_use, in_use, flux_use, flag_use))

    # Export CSV (python writer; no lock issues)
    export_table_to_csv(
        flux_table,
        CSV_OUT,
        preferred=["SLINK_ID", "Qmin_in", "Qmin_out", "Qmin_in_fix", "Qmin_out_fix", "Qmin_flux", "FLAG"]
    )

    # Optional: polygons + shapefile + flux raster
    if MAKE_SHP or MAKE_RASTER:
        print("\n8) OPTIONAL: RasterToPolygon + Dissolve + Join...")
        if arcpy.Exists(poly_raw):  arcpy.management.Delete(poly_raw)
        if arcpy.Exists(poly_diss): arcpy.management.Delete(poly_diss)
        if arcpy.Exists(poly_join): arcpy.management.Delete(poly_join)

        arcpy.conversion.RasterToPolygon(SUBCATCH, poly_raw, "NO_SIMPLIFY", "VALUE")
        # gridcode -> SLINK_ID
        flds = [f.name for f in arcpy.ListFields(poly_raw)]
        gc = [f for f in flds if f.lower() == "gridcode"]
        if gc:
            arcpy.management.AlterField(poly_raw, gc[0], "SLINK_ID", "SLINK_ID")

        arcpy.management.Dissolve(poly_raw, poly_diss, ["SLINK_ID"], multi_part="MULTI_PART")
        arcpy.management.CopyFeatures(poly_diss, poly_join)
        arcpy.management.JoinField(poly_join, "SLINK_ID", flux_table, "SLINK_ID",
                                   ["Qmin_in_fix", "Qmin_out_fix", "Qmin_flux", "FLAG"])

    if MAKE_SHP:
        print("\n9) OPTIONAL: Export shapefile...")
        if arcpy.Exists(SHP_OUT): arcpy.management.Delete(SHP_OUT)
        arcpy.conversion.FeatureClassToFeatureClass(poly_join, OUT_DIR, os.path.basename(SHP_OUT))
        print("SHP:", SHP_OUT)

    if MAKE_RASTER:
        print("\n10) OPTIONAL: PolygonToRaster (Qmin_flux) ...")
        # enforce exact template alignment
        arcpy.env.snapRaster = SLINK
        arcpy.env.cellSize   = SLINK
        arcpy.env.extent     = SLINK

        if os.path.exists(FLUX_RAS_OUT):
            try:
                os.remove(FLUX_RAS_OUT)
            except Exception:
                pass

        arcpy.conversion.PolygonToRaster(
            in_features=poly_join,
            value_field="Qmin_flux",
            out_rasterdataset=FLUX_RAS_OUT,
            cell_assignment="CELL_CENTER",
            priority_field="NONE",
            cellsize=SLINK
        )
        print("FLUX RASTER:", FLUX_RAS_OUT)

    print("\nDONE.")

if __name__ == "__main__":
    main()
