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
    target_ds = gdal.GetDriverByName('MEM').Create('', x_res, y_res, gdal.GDT_Byte)
    target_ds.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
    geotransform= target_ds.GetGeoTransform()
    band = target_ds.GetRasterBand(1)
    band.SetNoDataValue(255)
    
    # Rasterize
    gdal.RasterizeLayer(target_ds, [1], source_layer, burn_values=[1])
    
    # Read as array
    array = band.ReadAsArray()
    
    # Convert array to point coordinates
    count = 0
    roadList = np.where(array == 1)
    multipoint = ogr.Geometry(ogr.wkbMultiPoint)
    for indexY in roadList[0]:
        indexX = roadList[1][count]
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
    
    target_ds = None
    # Remove temporary files
    #while os.path.exists(working_directory + os.path.sep + 'temp.tif'):
    #    os.remove(working_directory + os.path.sep + 'temp.tif')

    spatialRef.MorphToESRI()
    with open(to_fname[:-4]+'.prj', 'w') as f:
        f.write(spatialRef.ExportToWkt())

    os.chdir(cwd)

            
if __name__ == "__main__":

    polygon_fn = r"C:\Workspace\part0075\MDB modelling\testbox\data_build\GW_model_area_model_buffer_20000.shp"
    direc = r"C:\Workspace\part0075\MDB modelling\testbox\data_build"    
    to_fname = "pilot_points.shp"     
    # Define pixel_size which equals distance betweens points
    pixel_size = 10000
    
    # Open the data source and read in the extent
    source_ds = ogr.Open(polygon_fn)
    poly2points(source_ds, to_fname, density=pixel_size, working_directory=direc)    
    source_ds = None
    
    from geopandas import read_file
    import matplotlib.pyplot as plt
    
    gw_model_poly = read_file(os.path.join(direc, "GW_model_area_model.shp"))
    pilot_points = read_file(os.path.join(direc, "pilot_points.shp"))
    
    fig = plt.figure(figsize=(8,10))
    ax = fig.add_subplot(1, 1, 1, aspect='equal')
    
    gw_model_poly.plot(ax=ax, alpha=0.2)
    
    pilot_points_xy = pilot_points['geometry'][0].wkt[12:-1].split(",")
    pilot_points_xy = [(float(x.strip(" ").split(" ")[0]), float(x.strip(" ").split(" ")[1])) for x in pilot_points_xy]
    #pilot_points.plot(markersize=5, ax=ax)
    plt.scatter([x[0] for x in pilot_points_xy], [x[1] for x in pilot_points_xy])
