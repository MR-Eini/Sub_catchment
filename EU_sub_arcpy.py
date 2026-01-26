#!/usr/bin/env python3
import os
import math
import time
import arcpy
from arcpy.sa import Fill, FlowDirection, FlowAccumulation, Con, StreamLink, Watershed, Raster

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True

# -----------------------------
# INPUTS
# -----------------------------
DEM_IN = r"E:\pythonProject\Europe\PREP_90M\dem_burned.tif"  # use burned DEM if you have it
# DEM_IN = r"E:\pythonProject\Europe\PREP_90M\dem_90m_3035.tif"        # or original 90m

OUT_DIR = r"E:\pythonProject\Europe\OUT_ARCGIS_25M"
SCRATCH_DIR = os.path.join(OUT_DIR, "_scratch")
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(SCRATCH_DIR, exist_ok=True)

# Stream definition threshold (km²)
THRESH_KM2 = 25.0

# Outputs (use .tif so rasterio/other tools can read later)
DEM_FILLED = os.path.join(OUT_DIR, "dem_filled.tif")
FDIR       = os.path.join(OUT_DIR, "flowdir_d8.tif")
FACC       = os.path.join(OUT_DIR, "flowacc_cells.tif")
UAA_KM2    = os.path.join(OUT_DIR, "upstream_area_km2.tif")
STREAM     = os.path.join(OUT_DIR, "stream_1000km2.tif")
SLINK      = os.path.join(OUT_DIR, "streamlink.tif")
CATCH      = os.path.join(OUT_DIR, "subcatch_streamlink.tif")
CATCH_POLY = os.path.join(OUT_DIR, "subcatch_streamlink.gpkg")  # optional

# -----------------------------
# ENVIRONMENT
# -----------------------------
arcpy.env.workspace = OUT_DIR
arcpy.env.scratchWorkspace = SCRATCH_DIR
arcpy.env.snapRaster = DEM_IN
arcpy.env.cellSize = DEM_IN
arcpy.env.extent = DEM_IN

# optional: compression (may or may not apply to all writes depending on tool)
try:
    arcpy.env.compression = "LZW"
except Exception:
    pass

def step(msg):
    print("\n" + msg)

def timed(label, fn):
    t0 = time.perf_counter()
    out = fn()
    dt = (time.perf_counter() - t0) / 60.0
    print(f"{label} DONE in {dt:.2f} min")
    return out

def main():
    if not os.path.exists(DEM_IN):
        raise FileNotFoundError(f"DEM not found: {DEM_IN}")

    desc = arcpy.Describe(DEM_IN)
    print("DEM:", DEM_IN)
    print("CRS:", desc.spatialReference.name)
    print("Extent:", desc.extent.XMin, desc.extent.YMin, desc.extent.XMax, desc.extent.YMax)

    dem_r = Raster(DEM_IN)
    resx = float(dem_r.meanCellWidth)
    resy = float(dem_r.meanCellHeight)
    print("Cellsize:", resx, resy)

    # EPSG:3035 => meters => constant cell area
    cell_km2 = (abs(resx) * abs(resy)) / 1e6
    thresh_cells = int(math.ceil(THRESH_KM2 / cell_km2))
    print(f"Cell area (km²): {cell_km2:.6f}")
    print(f"Stream threshold: {THRESH_KM2} km²  ~= {thresh_cells} cells")

    t_all = time.perf_counter()

    # 1) Fill sinks (ArcGIS)
    step("1) Fill depressions (ArcGIS Fill)...")
    if not os.path.exists(DEM_FILLED):
        def _fill():
            filled = Fill(dem_r)
            filled.save(DEM_FILLED)
        timed("Fill", _fill)
    else:
        print("Using existing:", DEM_FILLED)

    # 2) Flow direction D8
    step("2) Flow direction (D8)...")
    if not os.path.exists(FDIR):
        def _fdir():
            fd = FlowDirection(Raster(DEM_FILLED), flow_direction_type="D8")
            fd.save(FDIR)
        timed("FlowDirection", _fdir)
    else:
        print("Using existing:", FDIR)

    # 3) Flow accumulation (cells)
    step("3) Flow accumulation (cells)...")
    if not os.path.exists(FACC):
        def _facc():
            # FLOAT is safer than INTEGER for very large accumulations
            fa = FlowAccumulation(Raster(FDIR), data_type="FLOAT")
            fa.save(FACC)
        timed("FlowAccumulation", _facc)
    else:
        print("Using existing:", FACC)

    # 4) Upstream area in km²
    step("4) Convert flowacc cells -> upstream area (km²)...")
    if not os.path.exists(UAA_KM2):
        def _uaa():
            fa = Raster(FACC)
            uaa = fa * cell_km2
            uaa.save(UAA_KM2)
        timed("UpstreamArea km2", _uaa)
    else:
        print("Using existing:", UAA_KM2)

    # 5) Stream raster by km² threshold
    step("5) Stream raster (UAA_km2 >= 1000)...")
    if not os.path.exists(STREAM):
        def _stream():
            uaa = Raster(UAA_KM2)
            s = Con(uaa >= THRESH_KM2, 1)
            s.save(STREAM)
        timed("Stream thresholding", _stream)
    else:
        print("Using existing:", STREAM)

    # 6) StreamLink
    step("6) StreamLink...")
    if not os.path.exists(SLINK):
        def _slink():
            sl = StreamLink(Raster(STREAM), Raster(FDIR))
            sl.save(SLINK)
        timed("StreamLink", _slink)
    else:
        print("Using existing:", SLINK)

    # 7) Watershed per stream link (this yields subcatchments between junctions)
    step("7) Watershed (by StreamLink IDs) -> subcatchments...")
    if not os.path.exists(CATCH):
        def _catch():
            c = Watershed(Raster(FDIR), Raster(SLINK))
            c.save(CATCH)
        timed("Watershed", _catch)
    else:
        print("Using existing:", CATCH)

    # 8) Optional polygons (can be huge)
    # step("8) RasterToPolygon (optional)...")
    # if not arcpy.Exists(CATCH_POLY):
    #     def _poly():
    #         arcpy.conversion.RasterToPolygon(
    #             in_raster=CATCH,
    #             out_polygon_features=CATCH_POLY,
    #             simplify="NO_SIMPLIFY",
    #             raster_field="Value"
    #         )
    #     timed("RasterToPolygon", _poly)

    dt_all = (time.perf_counter() - t_all) / 60.0
    print(f"\nTOTAL time: {dt_all:.2f} min")
    print("Outputs:")
    print("  Filled DEM:", DEM_FILLED)
    print("  FlowDir:", FDIR)
    print("  FlowAcc:", FACC)
    print("  UpArea km²:", UAA_KM2)
    print("  Stream:", STREAM)
    print("  StreamLink:", SLINK)
    print("  Subcatch:", CATCH)

if __name__ == "__main__":
    main()
