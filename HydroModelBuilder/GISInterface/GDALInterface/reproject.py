# Code derived from:
# https://jgomezdans.github.io/gdal_notes/reprojection.html
import os
import subprocess

import numpy as np
from osgeo import gdal, gdalconst, ogr, osr


"""
This contains functions for reprojecting:
1. Rasters from one projected coordinate system to another, and allowing options:
    a. create_copy
    b. path_dest
    c. pixel_size
    d. new bounding box

2. Layer files, allowing the following options:
    a. create_copy
    b. path_dest
"""


def gdal_error_handler(err_class, err_num, err_msg):
    """

    :param err_class:
        :param err_num:
    :param err_msg:
    """
    errtype = {
        gdal.CE_None: 'None',
        gdal.CE_Debug: 'Debug',
        gdal.CE_Warning: 'Warning',
        gdal.CE_Failure: 'Failure',
        gdal.CE_Fatal: 'Fatal'
    }
    err_msg = err_msg.replace('\n', ' ')
    err_class = errtype.get(err_class, 'None')
    print 'Error Number: %s' % (err_num)
    print 'Error Type: %s' % (err_class)
    print 'Error Message: %s' % (err_msg)


def coordinate_transform(pointX, pointY, coordTransform):
    """A function to transform a point on the horizontal plane using
    CoordinateTransformation object from GDAL osr.

    :param pointX:
    :param pointY:
    :param coordTransform:
    """

    # create a geometry from coordinates
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(pointX, pointY)
    # transform point
    point.Transform(coordTransform)

    return point.GetX(), point.GetY()


def reproject_raster(ds,  # raster dataset
                     pixel_spacing=None,
                     src_cs=None,
                     dst_cs=None,
                     reproj_method=gdal.GRA_Bilinear,
                     create_copy=False,
                     copy_dest=None,
                     raster_driver="GTiff",
                     set_bounds=None):
    """A function to reproject and resample a GDAL dataset from within
    Python. The idea here is to reproject from one system to another, as well
    as to change the pixel size. The procedure is slightly long-winded, but
    goes like this:

    1. Set up the two Spatial Reference systems.
    2. Open the original dataset, and get the geotransform
    3. Calculate bounds of new geotransform by projecting the UL corners
    4. Calculate the number of pixels with the new projection & spacing
    5. Create an in-memory raster dataset
    6. Perform the projection

    Creates a copy of the projected raster to file if create_copy == True

    Returns projected gdal raster object at the end

    :param ds:
    :param pixel_spacing:  (Default value = None)
    :param src_cs: Default value = None)
    :param dst_cs: Default value = None)
    :param reproj_method: Default value = gdal.GRA_Bilinear)
    :param create_copy: Default value = False)
    :param copy_dest: Default value = None)
    :param raster_driver: Default value = "GTiff")
    :param set_bounds: Default value = None)
    """
    # Define the projection for "from", and get from dataset if not provided
    if src_cs == None:
        # Check raster GCS
        prj = ds.GetProjection()
        src_cs = osr.SpatialReference(wkt=prj)

    tx = osr.CoordinateTransformation(src_cs, dst_cs)

    # Get the Geotransform vector
    geo_t = ds.GetGeoTransform()
    x_size = ds.RasterXSize  # Raster xsize
    y_size = ds.RasterYSize  # Raster ysize
    NO_DATA = ds.GetRasterBand(1).GetNoDataValue()

    # Work out the boundaries of the new dataset in the target projection
    (ulx, uly) = coordinate_transform(geo_t[0], geo_t[3], tx)
    (lrx, lry) = coordinate_transform(geo_t[0] + geo_t[1] * x_size,
                                      geo_t[3] + geo_t[5] * y_size, tx)
    # Now, we create an in-memory raster
    mem_drv = gdal.GetDriverByName('MEM')
    # The size of the raster is given the new projection and pixel spacing
    # Using the values we calculated above. Also, setting it to store one band
    # and to use Float32 data type.

    if pixel_spacing == None:  # use default spacing from input raster
        pixel_spacing = geo_t[1]

    if set_bounds != None:
        # Define extent from set_bounds if it exists
        (ulx, lrx, lry, uly) = set_bounds

    row = int((lrx - ulx) / pixel_spacing)
    col = int((uly - lry) / pixel_spacing)

    dest = mem_drv.Create('', row, col, 1, gdal.GDT_Float32)
    # Calculate the new geotransform
    new_geo = (ulx, pixel_spacing, geo_t[2],
               uly, geo_t[4], -pixel_spacing)

    # Set the geotransform
    dest.SetGeoTransform(new_geo)
    dest.SetProjection(dst_cs.ExportToWkt())
    dest.GetRasterBand(1).SetNoDataValue(NO_DATA)
    dest.GetRasterBand(1).WriteArray(np.ones((col, row)) * NO_DATA)
    dest.GetRasterBand(1).FlushCache()
    # Perform the projection/resampling
    gdal.ReprojectImage(ds, dest, src_cs.ExportToWkt(), dst_cs.ExportToWkt(), reproj_method)

    # reproj_ds = dest.ReadAsArray()
    # Let's save it as a GeoTIFF.
    if create_copy:
        driver = gdal.GetDriverByName(raster_driver)
        metadata = driver.GetMetadata()
        # Create a copy of the transformed raster for external use
        if metadata.has_key(gdal.DCAP_CREATECOPY) and metadata[gdal.DCAP_CREATECOPY] == 'YES':
            if os.path.exists(copy_dest):
                # driver.DeleteDataSource(copy_dest)
                os.remove(copy_dest)
            driver.CreateCopy(copy_dest, dest, 0)  # , [ 'TILED=YES', 'COMPRESS=PACKBITS' ] )
            #dst_ds = None
        else:
            print 'Driver %s does not support CreateCopy() method.' % raster_driver

    return dest


"""
see www.gdal.org/ogr2ogr.html
"""


def ogr2ogr(src_datasource_name, dst_datasource_name, layer_name_list=[], **kwargs):
    """

    :param src_datasource_name:
    :param dst_datasource_name:
    :param layer_name_list: (Default value = [])
    :param dst_datasource_name:
    :param **kwargs:
    """
    if os.path.exists(dst_datasource_name):
        os.remove(dst_datasource_name)

    extras = []
    for key in kwargs.keys():
        extras.append('-' + key)
        if isinstance(kwargs[key], (list, tuple)):
            for v in kwargs[key]:
                extras.append(str(v))
        else:
            extras.append(str(kwargs[key]))

    ret = subprocess.check_output([
        "ogr2ogr",
        "-skipfailures",
        "-overwrite",
    ] + extras + [
        dst_datasource_name,
        src_datasource_name
    ] + layer_name_list)

    return ret


def reproject_layer(lyr_src,
                    src_cs=None,
                    dst_cs=None,
                    create_copy=False,
                    copy_dest=None,
                    geom_type=None):
    """
    :param lyr_src:
    :param src_cs:  (Default value = None)
    :param dst_cs: (Default value = None)
    :param create_copy: (Default value = False)
    :param copy_dest: (Default value = None)
    :param geom_type: (Default value = None)
    """

    ogr2ogr(src_datasource_name=lyr_src.GetName(),
            dst_datasource_name=copy_dest,
            t_srs=dst_cs,
            )
    driver = ogr.GetDriverByName("ESRI Shapefile")
    ds = driver.Open(copy_dest)
    return ds


if __name__ == "__main__":

    gdal.PushErrorHandler(gdal_error_handler)
    gdal.UseExceptions()
    # Test for raster reprojection
    dataset = r"C:\Workspace\part0075\MDB modelling\Campaspe_model\GIS\GIS_preprocessed\Hydrogeological_Unit_Layers\qa_1t_bb"
    pixel_spacing = 100
    epsg_from = None
    epsg_to = 28355
    reproj_method = gdal.GRA_Bilinear
    create_copy = True
    copy_dest = r"C:\Workspace\part0075\MDB modelling\test.tif"
    raster_driver = "EHdr"  # "GTiff"

    ds = gdal.Open(dataset, gdalconst.GA_ReadOnly)

    Proj_CS = osr.SpatialReference()
    # 28355 is the code for gda94 mga zone 55;
    # http://spatialreference.org/ref/epsg/gda94-mga-zone-55/
    Proj_CS.ImportFromEPSG(epsg_to)

    #set_corner=[ulx, uly, 0, lrx, lry, 0]
    reproject_raster(ds,  # raster dataset
                     pixel_spacing=pixel_spacing,
                     src_cs=None,
                     dst_cs=Proj_CS,
                     reproj_method=gdal.GRA_NearestNeighbour,
                     create_copy=create_copy,
                     copy_dest=copy_dest,
                     raster_driver="GTiff",
                     set_bounds=None)

    # Test for layer reprojection
    driver = ogr.GetDriverByName("ESRI Shapefile")
    ds = driver.Open(
        r"C:\Workspace\part0075\MDB modelling\testbox\input_data\Waterways\River_Murray.shp", 0)
    # r"C:\Workspace\part0075\MDB modelling\Campaspe_model\GIS\GIS_preprocessed\Surface_Water\Streams\Campaspe_Riv.shp", 0)
    poly_obj = ds.GetLayer()
    if poly_obj == None:
        print 'Could not open '

    srs = poly_obj.GetSpatialRef()
    # print srs.ExportToWkt()
    copy_file = "test_model.shp"
    copy_path = r"C:\Workspace\part0075\\MDB modelling"
    copy_dest = os.path.join(copy_path, copy_file)

    ds = reproject_layer(ds,
                         src_cs=srs,
                         dst_cs=Proj_CS,
                         create_copy=True,
                         copy_dest=copy_dest,
                         geom_type=ogr.wkbMultiLineString)
