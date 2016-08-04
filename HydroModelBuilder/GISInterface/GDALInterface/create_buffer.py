# -*- coding: utf-8 -*-
"""
Taken from http://snipplr.com/view/66633/
"""
from osgeo import ogr
import sys

def create_buffer4poly(shpfile, buffile=None, buffer_distance=1000.0):

    shp = ogr.Open(shpfile)
    if shp == None:
        print 'Could not open: ' + shpfile
        sys.exit(1)

    if buffile == None:
        buffile = shpfile[0:-4] + "_buffer_" + str(buffer_distance) + ".shp"
     # End if
        
    drv = shp.GetDriver()
    drv.CopyDataSource(shp, buffile)
     
    shp = None
     
    buf = ogr.Open(buffile, 1)
    lyr = buf.GetLayer(0)
    count=lyr.GetFeatureCount()
    
    for i in range(count):
        feat = lyr.GetFeature(i)
        lyr.DeleteFeature(i)
        geom = feat.GetGeometryRef()
        feat.SetGeometry(geom.Buffer(buffer_distance))
        lyr.CreateFeature(feat)

    buf.ExecuteSQL('REPACK ' + buffile.split("\\")[-1][0:-4])
    #buf.Destroy()

    return buf

if __name__ == '__main__':
    shpfile = r"C:\Workspace\part0075\MDB modelling\integrated\Modules\Groundwater\model_files\GW_model_area_model.shp"
    buffile = r"C:\Workspace\part0075\MDB modelling\integrated\Modules\Groundwater\model_files\GW_model_area_model_buffered.shp"
    buffer_distance = 10000.0
    create_buffer4poly(shpfile, buffile=buffile, buffer_distance=buffer_distance)    