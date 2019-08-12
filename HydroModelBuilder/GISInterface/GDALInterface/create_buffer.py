# -*- coding: utf-8 -*-
"""
Taken from http://snipplr.com/view/66633/
"""
import os
import sys
import warnings

from osgeo import ogr


def create_buffer4poly(shpfile, buffile=None, buffer_distance=1000.0):
    """
    :param shpfile:
    :param buffile:  (Default value = None)
    :param buffer_distance: (Default value = 1000.0)
    """
    shp = ogr.Open(shpfile)
    if shp == None:
        print('Could not open: ' + shpfile)
        sys.exit(1)

    if buffile == None:
        fname = os.path.splitext(os.path.basename(shpfile))[0]
        buffile = fname + "_buffer_" + str(buffer_distance) + ".shp"
    # End if

    # OGR fails to handle unicode characters!
    if type(buffile) is str:
        warnings.warn("unicode characters encountered - OGR sometimes fails to handle unicode")
        buffile = str(buffile)

    drv = shp.GetDriver()
    drv.CopyDataSource(shp, buffile)

    shp = None

    buf = ogr.Open(buffile, 1)
    lyr = buf.GetLayer(0)
    count = lyr.GetFeatureCount()

    for i in range(count):
        feat = lyr.GetFeature(i)
        lyr.DeleteFeature(i)
        geom = feat.GetGeometryRef()
        feat.SetGeometry(geom.Buffer(buffer_distance))
        lyr.CreateFeature(feat)

    buf.ExecuteSQL("REPACK " + os.path.splitext(buffile)[0])
    buf.FlushCache()
    # buf.Destroy()

    return buf
# End create_buffer4poly()


if __name__ == '__main__':
    shpfile = r"C:\Workspace\part0075\MDB modelling\integrated\Modules\Groundwater\model_files\GW_model_area_model.shp"
    buffile = r"C:\Workspace\part0075\MDB modelling\integrated\Modules\Groundwater\model_files\GW_model_area_model_buffered.shp"
    buffer_distance = 10000.0
    create_buffer4poly(shpfile, buffile=buffile, buffer_distance=buffer_distance)
