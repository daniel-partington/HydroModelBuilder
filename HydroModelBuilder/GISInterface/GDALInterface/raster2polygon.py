"""
Taken from: https://pcjericks.github.io/py-gdalogr-cookbook/raster_layers.html
"""

import sys

from osgeo import gdal, ogr

# this allows GDAL to throw Python Exceptions
gdal.UseExceptions()


def raster2polygon(raster, band_num=1):
    """

    :param raster:
    :param band_num:  (Default value = 1)
    """
    src_ds = gdal.Open(raster)
    if src_ds is None:
        print('Unable to open %s' % raster)
        sys.exit(1)

    try:
        srcband = src_ds.GetRasterBand(band_num)
    except RuntimeError as e:
        # for example, try GetRasterBand(10)
        print('Band ( %i ) not found' % band_num)
        print(e)
        sys.exit(1)

    #
    #  create output datasource
    #
    dst_layername = raster + "_polygon"
    drv = ogr.GetDriverByName("ESRI Shapefile")
    dst_ds = drv.CreateDataSource(dst_layername + ".shp")
    dst_layer = dst_ds.CreateLayer(dst_layername, srs=None)

    gdal.Polygonize(srcband, None, dst_layer, -1, [], callback=None)


if __name__ == "__main__":

    raster = r"C:\Workspace\part0075\MDB modelling\testbox\data_build\qa_1t_clipped.bil"
    raster2polygon(raster)

    from geopandas import read_file

    gw_model_poly = read_file(raster + '_polygon.shp')
    gw_model_poly.plot(alpha=0.2)
