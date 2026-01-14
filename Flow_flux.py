#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EU subcatchment flux balance from SSQ (values only on reaches)

Goal:
  For each subcatchment / streamlink ID, compute:
      inflow_sum = sum(SSQ_out of all immediate upstream links draining into this link)
      outflow    = SSQ_out at this link's downstream outlet
      deltaQ     = outflow - inflow_sum   (default sign convention)

This supports multiple inlets naturally.

Inputs (ArcGIS outputs on the SAME GRID, 90 m, EPSG:3035 recommended):
  - flowdir_d8.tif                 (ArcGIS FlowDirection, D8 codes)
  - flowacc_cells.tif              (ArcGIS FlowAccumulation, in CELLS)
  - streamlink.tif                 (ArcGIS StreamLink IDs)
  - subcatch_streamlink.tif         (ArcGIS Watershed by StreamLink IDs)
  - SSQ.tif                         (discharge/runoff on reaches; NoData or 0 elsewhere)

Outputs:
  - outlets (points) in FileGDB
  - shifted outlets (points) with DOWN_ID (downstream link id)
  - deltaQ_by_link.csv
  - deltaQ_table (for reclass)
  - deltaQ_raster.tif (deltaQ value for each subcatchment cell)

Run with ArcGIS Pro python (arcgispro-py3) in PyCharm.
"""

import os
import time
import csv
from pathlib import Path

import arcpy
from arcpy.sa import (
    Raster, SetNull, IsNull, Con, ZonalStatistics,
    ExtractMultiValuesToPoints, ExtractValuesToPoints,
    ReclassByTable
)

# -----------------------------
# USER SETTINGS
# -----------------------------
BASE = r"E:\pythonProject\Europe\OUT_ARCGIS_90M"   # folder containing your ArcGIS outputs
SSQ_PATH = r"E:\pythonProject\Europe\DELTAQ_90M\SSQ_90m_3035.tif"     # SSQ raster (values on reaches)

FLOWDIR   = os.path.join(BASE, "flowdir_d8.tif")
FLOWACC   = os.path.join(BASE, "flowacc_cells.tif")
STREAMLINK= os.path.join(BASE, "streamlink.tif")
SUBCATCH  = os.path.join(BASE, "subcatch_streamlink.tif")

OUT_DIR = r"E:\pythonProject\Europe\DELTAQ_90M"
Path(OUT_DIR).mkdir(parents=True, exist_ok=True)

# If SSQ uses 0 on non-reach cells, treat <=0 as "no data"
SSQ_NONREACH_IS_ZERO = True
SSQ_MIN_VALID = 0.0   # valid if > SSQ_MIN_VALID when SSQ_NONREACH_IS_ZERO is True

# If SSQ at the true streamlink outlet cell is missing, use fallback:
# pick the farthest-downstream cell that HAS SSQ inside each streamlink.
USE_SSQ_FALLBACK_OUTLET = True

# delta sign convention:
#   "out_minus_in" => deltaQ = outflow - sum(inflows)   (default; net lateral gain inside catchment)
#   "in_minus_out" => deltaQ = sum(inflows) - outflow
DELTA_SIGN = "out_minus_in"

# Parallel processing if available (can speed up raster ops)
PARALLEL = "75%"

DELTAQ_SCALE_DESIRED = 1000  # keeps 0.001 precision (deltaQ * 1000 stored as LONG)

# -----------------------------
# ArcGIS environment
# -----------------------------
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True
arcpy.env.workspace = OUT_DIR
arcpy.env.scratchWorkspace = os.path.join(OUT_DIR, "_scratch")
Path(arcpy.env.scratchWorkspace).mkdir(parents=True, exist_ok=True)
arcpy.env.parallelProcessingFactor = PARALLEL

# Use a common snap/cellsize/extent (critical)
arcpy.env.snapRaster = STREAMLINK
arcpy.env.cellSize = STREAMLINK
arcpy.env.extent = STREAMLINK


# -----------------------------
# Helpers
# -----------------------------
def now():
    return time.perf_counter()

def minutes(dt):
    return dt / 60.0

def log(msg):
    print(msg, flush=True)

def assert_exists(p, label):
    if not (p and arcpy.Exists(p) or os.path.exists(p)):
        raise FileNotFoundError(f"{label} not found: {p}")

def describe_raster(path):
    d = arcpy.Describe(path)
    sr = d.spatialReference
    return {
        "path": path,
        "sr_name": sr.name if sr else None,
        "sr_factory": getattr(sr, "factoryCode", None) if sr else None,
        "meanCellWidth": d.meanCellWidth,
        "meanCellHeight": d.meanCellHeight,
        "extent": (d.extent.XMin, d.extent.YMin, d.extent.XMax, d.extent.YMax),
        "rows": d.height,
        "cols": d.width
    }

def check_alignment(template, other, name_other):
    a = describe_raster(template)
    b = describe_raster(other)
    issues = []
    if a["sr_factory"] != b["sr_factory"]:
        issues.append(f"CRS differs: template EPSG {a['sr_factory']} vs {name_other} EPSG {b['sr_factory']}")
    if (abs(a["meanCellWidth"] - b["meanCellWidth"]) > 1e-6) or (abs(a["meanCellHeight"] - b["meanCellHeight"]) > 1e-6):
        issues.append(f"Cellsize differs: template {a['meanCellWidth']},{a['meanCellHeight']} vs {name_other} {b['meanCellWidth']},{b['meanCellHeight']}")
    # extent can differ slightly; we warn only if very different
    if a["rows"] != b["rows"] or a["cols"] != b["cols"]:
        issues.append(f"Rows/Cols differ: template {a['rows']}x{a['cols']} vs {name_other} {b['rows']}x{b['cols']}")
    if issues:
        log("WARNING alignment check:")
        for s in issues:
            log("  - " + s)
        log("If SSQ is not aligned, you should ProjectRaster/Resample it to match STREAMLINK before proceeding.\n")

def make_gdb(out_dir, name="deltaQ.gdb"):
    gdb = os.path.join(out_dir, name)
    if not arcpy.Exists(gdb):
        arcpy.management.CreateFileGDB(out_dir, name)
    return gdb

def delete_if_exists(path):
    if arcpy.Exists(path):
        arcpy.management.Delete(path)
    elif os.path.exists(path):
        try:
            os.remove(path)
        except Exception:
            pass

def delete_identical_by_field(fc, field):
    # Keep the first feature for each field value
    arcpy.management.DeleteIdentical(fc, [field])

def add_field_if_missing(fc, field, ftype):
    existing = {f.name.upper() for f in arcpy.ListFields(fc)}
    if field.upper() not in existing:
        arcpy.management.AddField(fc, field, ftype)

def calc_link_downstep_dxdy(flowdir_code, cell):
    # ArcGIS D8 FlowDirection codes:
    # 1=E, 2=SE, 4=S, 8=SW, 16=W, 32=NW, 64=N, 128=NE
    step = {
        1:   ( cell, 0.0),
        2:   ( cell,-cell),
        4:   ( 0.0,-cell),
        8:   (-cell,-cell),
        16:  (-cell, 0.0),
        32:  (-cell, cell),
        64:  ( 0.0, cell),
        128: ( cell, cell),
    }
    return step.get(int(flowdir_code), (0.0, 0.0))


def raster_outlet_points_from_streamlink_outlet(
    gdb, streamlink_ras, flowacc_ras, flowdir_ras, ssq_ras,
    tag="primary"
):
    """
    Outlet = downstream-most cell of each streamlink (max flowacc on streamlink cells).
    Creates outlet points with:
      LINK_ID, SSQ_OUT, FDIR
    """
    out_fc = os.path.join(gdb, f"outlets_{tag}")
    delete_if_exists(out_fc)

    # flowacc only where streamlink exists
    facc_on_stream = SetNull(IsNull(Raster(streamlink_ras)), Raster(flowacc_ras))
    facc_max = ZonalStatistics(streamlink_ras, "Value", facc_on_stream, "MAXIMUM", "DATA")
    outlet_mask = (facc_on_stream == facc_max)

    # store link id at outlet cells
    outlet_id_ras = Con(outlet_mask, Raster(streamlink_ras))
    outlet_id_ras_path = os.path.join(OUT_DIR, f"outlet_id_{tag}.tif")
    delete_if_exists(outlet_id_ras_path)
    outlet_id_ras.save(outlet_id_ras_path)

    # raster to point -> grid_code is link id
    arcpy.conversion.RasterToPoint(outlet_id_ras_path, out_fc, "VALUE")

    # remove duplicates (ties) per link id (grid_code)
    delete_identical_by_field(out_fc, "grid_code")

    # add LINK_ID
    add_field_if_missing(out_fc, "LINK_ID", "LONG")
    arcpy.management.CalculateField(out_fc, "LINK_ID", "!grid_code!", "PYTHON3")

    # sample SSQ and FDIR at the outlet point
    ExtractMultiValuesToPoints(out_fc, [[ssq_ras, "SSQ_OUT"], [flowdir_ras, "FDIR"]], "NONE")

    return out_fc


def raster_outlet_points_from_farthest_ssq_cell(
    gdb, streamlink_ras, flowacc_ras, flowdir_ras, ssq_ras,
    tag="fallback"
):
    """
    Fallback outlet per streamlink:
      pick farthest downstream cell that HAS SSQ (and is on streamlink), based on max flowacc.
    """
    out_fc = os.path.join(gdb, f"outlets_{tag}")
    delete_if_exists(out_fc)

    ssq = Raster(ssq_ras)
    sl  = Raster(streamlink_ras)
    facc= Raster(flowacc_ras)

    # define "SSQ exists on stream"
    if SSQ_NONREACH_IS_ZERO:
        ssq_valid = (~IsNull(ssq)) & (ssq > SSQ_MIN_VALID)
    else:
        ssq_valid = ~IsNull(ssq)

    mask = ssq_valid & (~IsNull(sl))
    facc_on_ssq = SetNull(~mask, facc)

    facc_max = ZonalStatistics(streamlink_ras, "Value", facc_on_ssq, "MAXIMUM", "DATA")
    outlet_mask = (facc_on_ssq == facc_max)

    outlet_id_ras = Con(outlet_mask, sl)
    outlet_id_ras_path = os.path.join(OUT_DIR, f"outlet_id_{tag}.tif")
    delete_if_exists(outlet_id_ras_path)
    outlet_id_ras.save(outlet_id_ras_path)

    arcpy.conversion.RasterToPoint(outlet_id_ras_path, out_fc, "VALUE")
    delete_identical_by_field(out_fc, "grid_code")

    add_field_if_missing(out_fc, "LINK_ID", "LONG")
    arcpy.management.CalculateField(out_fc, "LINK_ID", "!grid_code!", "PYTHON3")

    ExtractMultiValuesToPoints(out_fc, [[ssq_ras, "SSQ_OUT"], [flowdir_ras, "FDIR"]], "NONE")

    return out_fc


def shift_points_one_cell_downstream(gdb, outlet_fc, cellsize, tag="shifted"):
    """
    Create shifted points by moving each point 1 cell downstream using its FDIR.
    """
    sr = arcpy.Describe(outlet_fc).spatialReference
    out_fc = os.path.join(gdb, f"outlets_{tag}")
    delete_if_exists(out_fc)

    arcpy.management.CreateFeatureclass(
        out_path=os.path.dirname(out_fc),
        out_name=os.path.basename(out_fc),
        geometry_type="POINT",
        spatial_reference=sr
    )
    add_field_if_missing(out_fc, "LINK_ID", "LONG")

    # we keep a stable join key
    add_field_if_missing(out_fc, "OID_SRC", "LONG")

    # read original and write shifted
    with arcpy.da.SearchCursor(outlet_fc, ["OID@", "LINK_ID", "FDIR", "SHAPE@XY"]) as sc, \
         arcpy.da.InsertCursor(out_fc, ["OID_SRC", "LINK_ID", "SHAPE@XY"]) as ic:
        for oid, link_id, fdir, (x, y) in sc:
            dx, dy = calc_link_downstep_dxdy(fdir, cellsize)
            ic.insertRow((oid, link_id, (x + dx, y + dy)))

    return out_fc


def sample_downstream_linkid(gdb, shifted_fc, streamlink_ras, tag="with_down"):
    """
    Sample STREAMLINK at shifted points -> DOWN_ID.
    """
    out_fc = os.path.join(gdb, f"{os.path.basename(shifted_fc)}_{tag}")
    delete_if_exists(out_fc)

    ExtractValuesToPoints(shifted_fc, streamlink_ras, out_fc, "NONE", "VALUE_ONLY")

    # ExtractValuesToPoints writes sampled value into field "RASTERVALU"
    add_field_if_missing(out_fc, "DOWN_ID", "LONG")

    # DOWN_ID = int(RASTERVALU) else 0
    expr = "0 if !RASTERVALU! is None else int(!RASTERVALU!)"
    arcpy.management.CalculateField(out_fc, "DOWN_ID", expr, "PYTHON3")

    return out_fc


def read_outlet_attributes(out_fc):
    """
    Return dict: link_id -> (q_out, down_id)
    """
    q_out = {}
    down_id = {}

    # Detect field presence
    fields = [f.name for f in arcpy.ListFields(out_fc)]
    need = ["LINK_ID", "SSQ_OUT", "DOWN_ID"]
    for n in need:
        if n not in fields:
            raise RuntimeError(f"Missing field {n} in {out_fc}. Available: {fields}")

    with arcpy.da.SearchCursor(out_fc, ["LINK_ID", "SSQ_OUT", "DOWN_ID"]) as cur:
        for lid, q, did in cur:
            if lid is None:
                continue
            lid = int(lid)
            # SSQ_OUT might be None if SSQ missing at that outlet
            q_out[lid] = None if q is None else float(q)
            down_id[lid] = 0 if did is None else int(did)

    return q_out, down_id


def compute_deltaQ(q_out_by_id, down_by_id, sign="out_minus_in"):
    """
    inflow_sum[k] = sum(q_out[i] for i with down_id==k)
    deltaQ[k] = outflow[k] - inflow_sum[k]  (or reversed)
    """
    # fill missing q_out with 0 (but keep track)
    missing = [k for k, v in q_out_by_id.items() if v is None]
    q = {k: (0.0 if v is None else float(v)) for k, v in q_out_by_id.items()}

    inflow_sum = {k: 0.0 for k in q.keys()}
    for up_id, qout in q.items():
        dn = down_by_id.get(up_id, 0)
        if dn in inflow_sum:
            inflow_sum[dn] += qout

    delta = {}
    for k in q.keys():
        if sign == "out_minus_in":
            delta[k] = q[k] - inflow_sum.get(k, 0.0)
        elif sign == "in_minus_out":
            delta[k] = inflow_sum.get(k, 0.0) - q[k]
        else:
            raise ValueError("sign must be 'out_minus_in' or 'in_minus_out'")

    return q, inflow_sum, delta, missing


def write_csv(path, q, inflow_sum, delta, down_by_id):
    ids = sorted(q.keys())
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["id", "down_id", "q_out", "q_in_sum", "deltaQ"])
        for k in ids:
            w.writerow([k, down_by_id.get(k, 0), q.get(k, 0.0), inflow_sum.get(k, 0.0), delta.get(k, 0.0)])


def _choose_scale(delta_by_id, desired=1000, max_int=2_000_000_000):
    vals = [abs(float(v)) for v in delta_by_id.values() if v is not None]
    if not vals:
        return 1
    max_abs = max(vals)
    if max_abs == 0:
        return 1
    cap = int(max_int // max_abs)  # max scale that won't overflow LONG
    return max(1, min(int(desired), cap))

def build_reclass_table(gdb, table_name, delta_by_id, desired_scale=1000):
    """
    ReclassByTable requires integer output field => store scaled deltaQ as LONG.
    Returns: (table_path, from_field, to_field, new_int_field, scale)
    """
    tbl = os.path.join(gdb, table_name)
    delete_if_exists(tbl)

    arcpy.management.CreateTable(gdb, table_name)

    from_f = "FROM_ID"
    to_f   = "TO_ID"
    new_f  = "NEW_I"   # integer

    arcpy.management.AddField(tbl, from_f, "LONG")
    arcpy.management.AddField(tbl, to_f,   "LONG")
    arcpy.management.AddField(tbl, new_f,  "LONG")

    scale = _choose_scale(delta_by_id, desired=desired_scale)

    with arcpy.da.InsertCursor(tbl, [from_f, to_f, new_f]) as ic:
        for k, v in delta_by_id.items():
            kk = int(k)
            vv = 0.0 if v is None else float(v)
            ic.insertRow((kk, kk, int(round(vv * scale))))

    return tbl, from_f, to_f, new_f, scale


def align_raster_to_template(ssq_in, template_ras, out_ras, resampling="NEAREST"):
    """
    Force ssq_in onto the exact grid of template_ras:
      - same CRS
      - same cellsize
      - same extent (rows/cols)
      - snapped to template grid

    resampling:
      - "NEAREST" recommended here because SSQ values exist only on reaches and you want to preserve them.
    """
    tpl = arcpy.Describe(template_ras)
    tpl_sr = tpl.spatialReference
    tpl_cell = tpl.meanCellWidth
    tpl_ext = tpl.extent
    rect = f"{tpl_ext.XMin} {tpl_ext.YMin} {tpl_ext.XMax} {tpl_ext.YMax}"

    # If already aligned (same sr + same cell + same rows/cols), return input
    a = describe_raster(template_ras)
    b = describe_raster(ssq_in)
    same_sr = (a["sr_factory"] == b["sr_factory"])
    same_cell = (abs(a["meanCellWidth"] - b["meanCellWidth"]) < 1e-6 and abs(a["meanCellHeight"] - b["meanCellHeight"]) < 1e-6)
    same_rc = (a["rows"] == b["rows"] and a["cols"] == b["cols"])
    if same_sr and same_cell and same_rc:
        log("SSQ already matches template rows/cols, cellsize, CRS. No alignment needed.")
        return ssq_in

    tmp1 = os.path.join(arcpy.env.scratchWorkspace, "ssq_tmp1.tif")
    tmp2 = os.path.join(arcpy.env.scratchWorkspace, "ssq_tmp2.tif")
    delete_if_exists(tmp1)
    delete_if_exists(tmp2)
    delete_if_exists(out_ras)

    # Ensure snap/cellsize are set to template during processing
    old_snap = arcpy.env.snapRaster
    old_cell = arcpy.env.cellSize
    old_ext  = arcpy.env.extent

    arcpy.env.snapRaster = template_ras
    arcpy.env.cellSize = template_ras
    arcpy.env.extent = template_ras  # many tools respect it; we still Clip explicitly

    try:
        # 1) Project if CRS differs
        ssq_desc = arcpy.Describe(ssq_in)
        ssq_sr = ssq_desc.spatialReference
        if (not ssq_sr) or (ssq_sr.factoryCode != tpl_sr.factoryCode):
            log(f"Align SSQ: projecting to EPSG:{tpl_sr.factoryCode} ...")
            arcpy.management.ProjectRaster(
                in_raster=ssq_in,
                out_raster=tmp1,
                out_coor_system=tpl_sr,
                resampling_type=resampling,
                cell_size=tpl_cell
            )
            src_for_resample = tmp1
        else:
            src_for_resample = ssq_in

        # 2) Resample if cellsize differs
        src_desc = arcpy.Describe(src_for_resample)
        if abs(src_desc.meanCellWidth - tpl_cell) > 1e-6 or abs(src_desc.meanCellHeight - tpl_cell) > 1e-6:
            log(f"Align SSQ: resampling to {tpl_cell} ...")
            arcpy.management.Resample(
                in_raster=src_for_resample,
                out_raster=tmp2,
                cell_size=tpl_cell,
                resampling_type=resampling
            )
            src_for_clip = tmp2
        else:
            src_for_clip = src_for_resample

        # 3) Clip to template extent to force identical rows/cols
        log("Align SSQ: clipping to template extent (forces rows/cols match) ...")
        arcpy.management.Clip(
            in_raster=src_for_clip,
            rectangle=rect,
            out_raster=out_ras,
            in_template_dataset=None,
            nodata_value=None,
            clipping_geometry="NONE",
            maintain_clipping_extent="MAINTAIN_EXTENT"
        )

    finally:
        # restore env
        arcpy.env.snapRaster = old_snap
        arcpy.env.cellSize = old_cell
        arcpy.env.extent = old_ext

    # Final sanity check: must match rows/cols now
    check_alignment(template_ras, out_ras, "SSQ_ALIGNED")
    a2 = describe_raster(template_ras)
    b2 = describe_raster(out_ras)
    if a2["rows"] != b2["rows"] or a2["cols"] != b2["cols"]:
        raise RuntimeError(
            f"SSQ alignment failed: template {a2['rows']}x{a2['cols']} vs aligned {b2['rows']}x{b2['cols']}. "
            f"Check extents/CRS."
        )

    log("SSQ aligned raster: " + out_ras)
    return out_ras

# -----------------------------
# MAIN
# -----------------------------
def main():
    t_all = now()

    # 0) checks
    log("0) Checking inputs...")
    for p, lab in [
        (FLOWDIR, "FLOWDIR"),
        (FLOWACC, "FLOWACC"),
        (STREAMLINK, "STREAMLINK"),
        (SUBCATCH, "SUBCATCH"),
        (SSQ_PATH, "SSQ"),
    ]:
        assert_exists(p, lab)

    check_alignment(STREAMLINK, SSQ_PATH, "SSQ")

    # Force SSQ onto the STREAMLINK grid (recommended)
    ssq_aligned_path = os.path.join(OUT_DIR, "SSQ_aligned.tif")
    SSQ_USED = align_raster_to_template(SSQ_PATH, STREAMLINK, ssq_aligned_path, resampling="NEAREST")

    d = describe_raster(STREAMLINK)
    cell = float(d["meanCellWidth"])
    log(f"Template grid: cell={cell} m, rows/cols={d['rows']} x {d['cols']}, EPSG={d['sr_factory']} ({d['sr_name']})")

    # outputs workspace
    gdb = make_gdb(OUT_DIR, "deltaQ.gdb")
    log("Output GDB: " + gdb)

    # 1) outlet points (primary = true streamlink outlet)
    t1 = now()
    log("1) Building outlet points (primary: max flowacc on streamlink)...")
    outlets_primary = raster_outlet_points_from_streamlink_outlet(
        gdb, STREAMLINK, FLOWACC, FLOWDIR, SSQ_USED, tag="primary"
    )
    log(f"   outlets_primary: {outlets_primary}")
    log(f"   Done in {minutes(now()-t1):.2f} min")

    # 2) optional fallback outlet points (farthest downstream cell with SSQ)
    outlets_fallback = None
    if USE_SSQ_FALLBACK_OUTLET:
        t2 = now()
        log("2) Building outlet points (fallback: farthest downstream cell WITH SSQ on each streamlink)...")
        outlets_fallback = raster_outlet_points_from_farthest_ssq_cell(
            gdb, STREAMLINK, FLOWACC, FLOWDIR, SSQ_USED, tag="fallback"
        )
        log(f"   outlets_fallback: {outlets_fallback}")
        log(f"   Done in {minutes(now()-t2):.2f} min")

    # 3) shift primary points one cell downstream and sample downstream streamlink id
    t3 = now()
    log("3) Shifting outlet points 1 cell downstream and sampling DOWN_ID...")
    shifted = shift_points_one_cell_downstream(gdb, outlets_primary, cellsize=cell, tag="shifted")
    shifted_down = sample_downstream_linkid(gdb, shifted, STREAMLINK, tag="down")
    log(f"   shifted_down: {shifted_down}")
    log(f"   Done in {minutes(now()-t3):.2f} min")

    # 4) bring SSQ_OUT + DOWN_ID together in one FC (join by LINK_ID)
    #    Easiest: update shifted_down with SSQ_OUT from outlets_primary (by LINK_ID)
    t4 = now()
    log("4) Joining SSQ_OUT onto shifted_down (by LINK_ID)...")

    # Build SSQ map from outlets_primary
    ssq_map_primary = {}
    with arcpy.da.SearchCursor(outlets_primary, ["LINK_ID", "SSQ_OUT"]) as cur:
        for lid, q in cur:
            if lid is None:
                continue
            ssq_map_primary[int(lid)] = None if q is None else float(q)

    # Fallback SSQ if enabled
    ssq_map_fallback = {}
    if outlets_fallback:
        with arcpy.da.SearchCursor(outlets_fallback, ["LINK_ID", "SSQ_OUT"]) as cur:
            for lid, q in cur:
                if lid is None:
                    continue
                ssq_map_fallback[int(lid)] = None if q is None else float(q)

    add_field_if_missing(shifted_down, "SSQ_OUT", "DOUBLE")

    with arcpy.da.UpdateCursor(shifted_down, ["LINK_ID", "SSQ_OUT"]) as cur:
        for lid, q in cur:
            if lid is None:
                continue
            lid = int(lid)
            qq = ssq_map_primary.get(lid, None)
            if (qq is None) and USE_SSQ_FALLBACK_OUTLET:
                qq = ssq_map_fallback.get(lid, None)
            cur.updateRow([lid, qq])

    log(f"   Done in {minutes(now()-t4):.2f} min")

    # 5) compute deltaQ
    t5 = now()
    log("5) Computing deltaQ per link/subcatchment...")

    # Build dicts from shifted_down
    q_out = {}
    down_id = {}
    missing_q = 0

    with arcpy.da.SearchCursor(shifted_down, ["LINK_ID", "SSQ_OUT", "DOWN_ID"]) as cur:
        for lid, q, did in cur:
            if lid is None:
                continue
            lid = int(lid)
            q_out[lid] = None if q is None else float(q)
            if q is None:
                missing_q += 1
            down_id[lid] = 0 if did is None else int(did)

    q, inflow_sum, delta, missing = compute_deltaQ(q_out, down_id, sign=DELTA_SIGN)

    log(f"   Links found: {len(q)}")
    log(f"   Missing SSQ_OUT (after fallback if enabled): {len(missing)}")
    log(f"   Done in {minutes(now()-t5):.2f} min")

    # 6) write CSV
    t6 = now()
    csv_out = os.path.join(OUT_DIR, "deltaQ_by_link.csv")
    log("6) Writing CSV: " + csv_out)
    write_csv(csv_out, q, inflow_sum, delta, down_id)
    log(f"   Done in {minutes(now()-t6):.2f} min")

    # 7) build deltaQ raster by reclassing SUBCATCH (IDs) using ReclassByTable
    t7 = now()
    log("7) Building deltaQ raster via ReclassByTable (SUBCATCH IDs -> deltaQ)...")
    tbl, from_f, to_f, new_f, scale = build_reclass_table(gdb, "deltaQ_table", delta,
                                                          desired_scale=DELTAQ_SCALE_DESIRED)

    delta_int = ReclassByTable(SUBCATCH, tbl, from_f, to_f, new_f, "NODATA")

    # convert scaled integers back to float deltaQ
    delta_ras = arcpy.sa.Float(delta_int) / float(scale)

    delta_ras_path = os.path.join(OUT_DIR, "deltaQ_subcatch.tif")
    delete_if_exists(delta_ras_path)
    delta_ras.save(delta_ras_path)

    log("   deltaQ raster: " + delta_ras_path)
    log(f"   Done in {minutes(now() - t7):.2f} min")

    # Done
    log(f"TOTAL time: {minutes(now()-t_all):.2f} min")
    log("Outputs:")
    log("  Points (GDB): " + gdb)
    log("  Shifted points with DOWN_ID + SSQ_OUT: " + shifted_down)
    log("  CSV: " + csv_out)
    log("  deltaQ raster: " + delta_ras_path)


if __name__ == "__main__":
    main()
