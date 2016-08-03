# -*- coding: utf-8 -*-
"""
Created on Thu Mar 03 13:51:01 2016

@author: part0075
"""

from osgeo import gdal, gdalconst, osr
import numpy as np


def create_basement_bottom(hu_raster_path, surface_raster_file, basement_top_raster_file, basement_bot_raster_file, output_path, raster_driver='GTiff'):
    # Read in data from surface to array
    
    surf = gdal.Open(hu_raster_path + surface_raster_file, gdalconst.GA_ReadOnly)
    basetop = gdal.Open(hu_raster_path + basement_top_raster_file, gdalconst.GA_ReadOnly)
    
    surf_array = surf.GetRasterBand(1).ReadAsArray()
    surf_NODATA = surf.GetRasterBand(1).GetNoDataValue
    basetop_array = basetop.GetRasterBand(1).ReadAsArray()
    basetop_NODATA = basetop.GetRasterBand(1).GetNoDataValue()
    basetop_geotransform = basetop.GetGeoTransform()
    basetop_proj = basetop.GetProjection()
    src_cs=osr.SpatialReference(wkt=basetop_proj)
    
    def thickness_rule(top, bottom):
        # Applying the same rule as in Hocking
        # i.e. 100m thick or 10m thick if >250m and linearly vary in between    
        min_thickness = 10.0    
        max_thickness = 100.0
        upper_bound = 250.0
        lower_bound = 0.0    
        
        diff = top - bottom
        if diff < 0:
            print 'Error, top is below bottom'
        if diff >= upper_bound:
            thickness = min_thickness
        elif diff == lower_bound:
            thickness = max_thickness
        elif diff > lower_bound and diff < upper_bound:
            thickness = np.interp(diff, [lower_bound, upper_bound], [min_thickness, max_thickness])
    
        return bottom-thickness        
    
    vecfunc = np.vectorize(thickness_rule, otypes=[np.float64])
    
    basebot_array = vecfunc(surf_array, basetop_array)
    # Create mask of data with NODATA
    mask = np.ma.masked_array(basetop_array, basetop_array==basetop_NODATA)
    basebot_array[mask.mask] = basetop_NODATA
    
    # Create basement raster with new array

    def array2raster(newRasterfn, geoTransform, NODATA, src_cs, array, raster_driver):
    
        cols = array.shape[1]
        rows = array.shape[0]
    
        driver = gdal.GetDriverByName(raster_driver)
        outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Float64)
        outRaster.SetGeoTransform(geoTransform)
        outband = outRaster.GetRasterBand(1)
        outband.WriteArray(array)
        outband.SetNoDataValue(array[0][0])
        outRaster.SetProjection(src_cs.ExportToWkt())
        outband.FlushCache()        

    array2raster(output_path+basement_bot_raster_file+'.tif', basetop_geotransform, basetop_NODATA, src_cs, basebot_array, raster_driver)    

if __name__ == "__main__":
    hu_raster_path = r"C:\Workspace\part0075\MDB modelling\VAF_v2.0_ESRI_GRID\ESRI_GRID\Preprocessed_data\Hydrogeological_Unit_Layers\\" #r"C:\Workspace\part0075\MDB modelling\Campaspe_model\GIS\GIS_preprocessed\Hydrogeological_Unit_Layers\\"
    output_path = hu_raster_path  
    surface_raster_file = "sur_1t_mb"
    basement_top_raster_file = "bse_1t_mb"
    basement_bot_raster_file = "bse_2b_mb"
    raster_driver = 'GTiff'
    create_basement_bottom(hu_raster_path, surface_raster_file, basement_top_raster_file, basement_bot_raster_file, output_path, raster_driver=raster_driver)