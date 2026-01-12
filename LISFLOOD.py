#!/usr/bin/env python3
"""
Delineate LISFLOOD-consistent subcatchments from LDD + Upstream Area maps.

Requirements:
    pip install rasterio geopandas shapely whitebox
"""

import rasterio
import numpy as np
import geopandas as gpd
from shapely.geometry import Point
import whitebox

# --------------------------------------------------------------------
#                 USER-DEFINED INPUTS (REPLACE THESE)
# --------------------------------------------------------------------
LDD_PATH = "/path/to/ldd_map.tif"          # LISFLOOD LDD raster
UAA_PATH = "/path/to/upstream_area.tif"    # Upstream area raster (m²)
THRESH_KM2 = 1000                          # threshold for river extraction
OUTLETS_SHP = "/path/to/outlets.shp"       # output outlets
SUBCATCH_RASTER = "/path/to/subcatch.tif"  # output subcatchment raster
SUBCATCH_SHP = "/path/to/subcatch.shp"     # output shapefile
WBT_BINARY = "/path/to/whitebox_tools"     # path to WhiteboxTools executable
# --------------------------------------------------------------------

wbt = whitebox.WhiteboxTools()
wbt.set_whitebox_dir(WBT_BINARY)

# -----------------------------------------------------------
# 1. Load UAA and LDD maps
# -----------------------------------------------------------
with rasterio.open(UAA_PATH) as ds:
    uaa = ds.read(1)
    profile = ds.profile
    transform = ds.transform
    nodata = ds.nodata
    height, width = uaa.shape

with rasterio.open(LDD_PATH) as ds_ldd:
    ldd = ds_ldd.read(1)

# Convert UAA from m² → km² (if needed)
uaa_km2 = uaa / 1e6

# -----------------------------------------------------------
# 2. Extract river network (threshold)
# -----------------------------------------------------------
river = (uaa_km2 >= THRESH_KM2)

# -----------------------------------------------------------
# 3. Function to determine upstream neighbors in LISFLOOD LDD
# -----------------------------------------------------------
# LISFLOOD uses PCRaster-style LDD:
# 1=E, 2=NE, 3=N, 4=NW, 5=W, 6=SW, 7=S, 8=SE

# Mapping of LDD codes to cell offsets (dy, dx)
ldd_dirs = {
    1: ( 0,  1),  # E
    2: (-1,  1),  # NE
    3: (-1,  0),  # N
    4: (-1, -1),  # NW
    5: ( 0, -1),  # W
    6: ( 1, -1),  # SW
    7: ( 1,  0),  # S
    8: ( 1,  1),  # SE
}

# Reverse mapping: for each cell, which directions could flow INTO it
reverse_neighbors = {}
for code, (dy, dx) in ldd_dirs.items():
    reverse_neighbors.setdefault((dy, dx), code)

# -----------------------------------------------------------
# 4. Identify headwaters & confluences
# -----------------------------------------------------------
outlet_points = []

for row in range(height):
    for col in range(width):

        if not river[row, col]:
            continue

        # Find how many neighbors drain into this cell
        upstream_count = 0
        for code, (dy, dx) in ldd_dirs.items():
            r2 = row + dy
            c2 = col + dx
            if 0 <= r2 < height and 0 <= c2 < width:
                # Does (r2, c2) flow into (row, col)?
                if ldd[r2, c2] == reverse_neighbors.get((-dy, -dx)):
                    upstream_count += 1

        # Headwater: upstream_count == 0
        # Confluence: upstream_count >= 2
        if upstream_count == 0 or upstream_count >= 2:
            # Convert row,col to map coordinates
            x, y = rasterio.transform.xy(transform, row, col)
            outlet_points.append(Point(x, y))

# -----------------------------------------------------------
# 5. Save outlet points as shapefile
# -----------------------------------------------------------
gdf = gpd.GeoDataFrame({"id": range(1, len(outlet_points)+1)},
                       geometry=outlet_points,
                       crs=profile["crs"])
gdf.to_file(OUTLETS_SHP)

print(f"Saved {len(outlet_points)} outlet points → {OUTLETS_SHP}")

# -----------------------------------------------------------
# 6. Run watershed delineation (WhiteboxTools)
# -----------------------------------------------------------
wbt.watershed(
    d8_pntr=LDD_PATH,
    pour_pts=OUTLETS_SHP,
    output=SUBCATCH_RASTER,
    esri_pntr=False  # LISFLOOD uses standard D8 pointer (PCRaster style)
)

print(f"Subcatchment raster written to: {SUBCATCH_RASTER}")

# -----------------------------------------------------------
# 7. Convert raster → polygons
# -----------------------------------------------------------
wbt.raster_to_vector_polygons(
    input=SUBCATCH_RASTER,
    output=SUBCATCH_SHP
)

print(f"Subcatchment polygons written to: {SUBCATCH_SHP}")

