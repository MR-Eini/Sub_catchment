import arcpy
snap = r"E:\pythonProject\Europe\OUT_ARCGIS_90M\streamlink.tif"
mmq  = r"E:\pythonProject\Europe\DELTAQ_90M_MMQ\MMQ_90m_3035.tif"

for p in [snap, mmq]:
    d = arcpy.Describe(p)
    csx = arcpy.management.GetRasterProperties(p, "CELLSIZEX").getOutput(0)
    csy = arcpy.management.GetRasterProperties(p, "CELLSIZEY").getOutput(0)
    print("\n", p)
    print("CRS:", d.spatialReference.factoryCode, d.spatialReference.name)
    print("Cell:", csx, csy)
    print("Extent:", d.extent.XMin, d.extent.YMin, d.extent.XMax, d.extent.YMax)
    print("Rows/Cols:", d.height, d.width)