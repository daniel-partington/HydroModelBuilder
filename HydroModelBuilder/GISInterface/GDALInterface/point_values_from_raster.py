# Based on snippet @ http://gis.stackexchange.com/questions/46893/how-do-i-get-the-pixel-value-of-a-gdal-raster-under-an-ogr-point-without-numpy

from osgeo import gdal,ogr

def get_point_values_from_raster(points, raster):
    gt=raster.GetGeoTransform()
    rb=raster.GetRasterBand(1)
    
    if isinstance(points, list):    
        point_vals = []
        for point in points:
            mx,my=point[0], point[1]  #coord in map units
    
            #Convert from map to pixel coordinates.
            #Only works for geotransforms with no rotation.
            px = int((mx - gt[0]) / gt[1]) #x pixel
            py = int((my - gt[3]) / gt[5]) #y pixel

            val=rb.ReadAsArray(px,py,1,1)
            point_vals += [val[0][0]] 

    else:
        lyr=points.GetLayer()
        point_vals = []
        for feat in lyr:
            geom = feat.GetGeometryRef()
            mx,my=geom.GetX(), geom.GetY()  #coord in map units
    
            #Convert from map to pixel coordinates.
            #Only works for geotransforms with no rotation.
            px = int((mx - gt[0]) / gt[1]) #x pixel
            py = int((my - gt[3]) / gt[5]) #y pixel

            val=rb.ReadAsArray(px,py,1,1)
            point_vals += [val[0][0]] 
        
        lyr = None
    
    return point_vals        
    
if __name__ == "__main__":
    
    points= 1
    rasters = 2

    src_filename = r'C:\Workspace\part0075\MDB modelling\integrated\Modules\Groundwater\model_files\utqa_1t_clipped.bil'
    shp_filename = r'C:\Workspace\part0075\MDB modelling\integrated\Modules\Groundwater\model_files\pumping wells_clipped.shp'

    src_ds=gdal.Open(src_filename) 
    ds=ogr.Open(shp_filename)
    
    points_values = get_point_values_from_raster(ds, src_ds)