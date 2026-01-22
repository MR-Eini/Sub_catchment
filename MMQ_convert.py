import os
import arcpy

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True

MMQ_IN = r"Y:\MapStatistic\hydrology\discharge\MMQ.tif"
SNAP   = r"E:\pythonProject\Europe\OUT_ARCGIS_90M\streamlink.tif"
OUT    = r"E:\pythonProject\Europe\DELTAQ_90M_MMQ_FIX"
os.makedirs(OUT, exist_ok=True)

sr_out = arcpy.Describe(SNAP).spatialReference
ext_snap = arcpy.Describe(SNAP).extent
cell = float(arcpy.management.GetRasterProperties(SNAP, "CELLSIZEX").getOutput(0))

# ---- Inspect original MMQ georeference (before defining anything) ----
d0 = arcpy.Describe(MMQ_IN)
cs0x = arcpy.management.GetRasterProperties(MMQ_IN, "CELLSIZEX").getOutput(0)
cs0y = arcpy.management.GetRasterProperties(MMQ_IN, "CELLSIZEY").getOutput(0)
print("ORIGINAL MMQ (as-is)")
print("  CRS:", d0.spatialReference.name, getattr(d0.spatialReference, "factoryCode", None))
print("  Cell:", cs0x, cs0y)
print("  Extent:", d0.extent.XMin, d0.extent.YMin, d0.extent.XMax, d0.extent.YMax)
print("  Rows/Cols:", d0.height, d0.width)

# ---- Define correct CRS on a copy (metadata only) ----
# Set this to the TRUE CRS of MMQ_IN. If the extent looks like degrees, EPSG:4326 is most likely.
ASSUMED_EPSG = 4326
MMQ_DEFINED = os.path.join(OUT, "MMQ_defined.tif")
arcpy.management.CopyRaster(MMQ_IN, MMQ_DEFINED)
arcpy.management.DefineProjection(MMQ_DEFINED, arcpy.SpatialReference(ASSUMED_EPSG))

sr_in = arcpy.SpatialReference(ASSUMED_EPSG)
transforms = arcpy.ListTransformations(sr_in, sr_out)
geo_transform = transforms[0] if transforms else None
print("Using transform:", geo_transform)

# ---- Project to 3035 with 90m + snap ----
arcpy.env.snapRaster = SNAP
arcpy.env.cellSize = cell
arcpy.env.outputCoordinateSystem = sr_out

MMQ_PRJ = os.path.join(OUT, "MMQ_prj_3035_90m.tif")
arcpy.management.ProjectRaster(
    in_raster=MMQ_DEFINED,
    out_raster=MMQ_PRJ,
    out_coor_system=sr_out,
    resampling_type="BILINEAR",
    cell_size=cell,
    geographic_transform=geo_transform
)

# ---- Clip to EXACT template extent (guarantee same rows/cols/extent) ----
rect = f"{ext_snap.XMin} {ext_snap.YMin} {ext_snap.XMax} {ext_snap.YMax}"
MMQ_FINAL = os.path.join(OUT, "MMQ_90m_3035_aligned.tif")

arcpy.management.Clip(
    in_raster=MMQ_PRJ,
    rectangle=rect,
    out_raster=MMQ_FINAL,
    in_template_dataset=SNAP,
    nodata_value="",
    clipping_geometry="NONE",
    maintain_clipping_extent="MAINTAIN_EXTENT"
)

print("Wrote:", MMQ_FINAL)

# ---- Report final properties ----
df = arcpy.Describe(MMQ_FINAL)
print("FINAL")
print("  CRS:", df.spatialReference.factoryCode, df.spatialReference.name)
print("  Cell:", arcpy.management.GetRasterProperties(MMQ_FINAL, "CELLSIZEX").getOutput(0),
              arcpy.management.GetRasterProperties(MMQ_FINAL, "CELLSIZEY").getOutput(0))
print("  Extent:", df.extent.XMin, df.extent.YMin, df.extent.XMax, df.extent.YMax)
print("  Rows/Cols:", df.height, df.width)
