import os
import numpy as np
from osgeo import ogr, gdal

def poly2points(polygon_obj, to_fname="points.shp", density=10, working_directory=os.getcwd()):
    """
    Function to produce points
    
    """

    cwd = os.getcwd()
    os.chdir(working_directory)
    
    source_layer = source_ds.GetLayer()
    spatialRef = source_layer.GetSpatialRef()
    x_min, x_max, y_min, y_max = source_layer.GetExtent()
    
    # Create the destination data source
    x_res = int((x_max - x_min) / pixel_size)
    y_res = int((y_max - y_min) / pixel_size)
    target_ds = gdal.GetDriverByName('GTiff').Create(working_directory + os.path.sep + 'temp.tif', x_res, y_res, gdal.GDT_Byte)
    target_ds.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
    band = target_ds.GetRasterBand(1)
    band.SetNoDataValue(255)
    
    # Rasterize
    gdal.RasterizeLayer(target_ds, [1], source_layer, burn_values=[1])
    
    # Read as array
    array = band.ReadAsArray()
    
    raster = gdal.Open(working_directory + os.path.sep + 'temp.tif')
    geotransform = raster.GetGeoTransform()
    
    # Convert array to point coordinates
    count = 0
    roadList = np.where(array == 1)
    multipoint = ogr.Geometry(ogr.wkbMultiPoint)
    for indexY in roadList[0]:
        indexX = roadList[1][count]
        geotransform = raster.GetGeoTransform()
        originX = geotransform[0]
        originY = geotransform[3]
        pixelWidth = geotransform[1]
        pixelHeight = geotransform[5]
        Xcoord = originX+pixelWidth*indexX
        Ycoord = originY+pixelHeight*indexY
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(Xcoord, Ycoord)
        multipoint.AddGeometry(point)
        count += 1
    
    # Write point coordinates to Shapefile
    shpDriver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(to_fname):
        shpDriver.DeleteDataSource(to_fname)
    outDataSource = shpDriver.CreateDataSource(to_fname)
    outLayer = outDataSource.CreateLayer(to_fname, geom_type=ogr.wkbMultiPoint)
    featureDefn = outLayer.GetLayerDefn()
    outFeature = ogr.Feature(featureDefn)
    outFeature.SetGeometry(multipoint)
    outLayer.CreateFeature(outFeature)
    
    raster =None
    target_ds = None
    # Remove temporary files
    while os.path.exists(working_directory + os.path.sep + 'temp.tif'):
        os.remove(working_directory + os.path.sep + 'temp.tif')

    spatialRef.MorphToESRI()
    with open(to_fname[:-4]+'.prj', 'w') as f:
        f.write(spatialRef.ExportToWkt())

    os.chdir(cwd)

            
if __name__ == "__main__":

    polygon_fn = r"C:\Workspace\part0075\MDB modelling\testbox\01_steady_state\GW_model_area_model.shp"
    direc = r"C:\Workspace\part0075\MDB modelling\testbox\01_steady_state\\"    
    to_fname = "points.shp"     
    # Define pixel_size which equals distance betweens points
    pixel_size = 5000
    
    # Open the data source and read in the extent
    source_ds = ogr.Open(polygon_fn)
    poly2points(source_ds, to_fname, density=pixel_size, working_directory=direc)    
    source_ds = None