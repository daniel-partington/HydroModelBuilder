import os, sys
import ogr
from math import ceil


def main(outputGridfn,xmin,xmax,ymin,ymax,gridHeight,gridWidth, spatialRef):
    """

    :param outputGridfn: param xmin:
    :param xmax: param ymin:
    :param ymax: param gridHeight:
    :param gridWidth: param spatialRef:
    :param xmin: 
    :param ymin: 
    :param gridHeight: 
    :param spatialRef: 

    """

    # convert sys.argv to float
    xmin = float(xmin)
    xmax = float(xmax)
    ymin = float(ymin)
    ymax = float(ymax)
    gridWidth = float(gridWidth)
    gridHeight = float(gridHeight)

    # get rows
    rows = ceil((ymax-ymin)/gridHeight)
    # get columns
    cols = ceil((xmax-xmin)/gridWidth)

    # start grid cell envelope
    ringXleftOrigin = xmin
    ringXrightOrigin = xmin + gridWidth
    ringYtopOrigin = ymax
    ringYbottomOrigin = ymax-gridHeight

    # create output file
    outDriver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(outputGridfn):
        os.remove(outputGridfn)
    outDataSource = outDriver.CreateDataSource(outputGridfn)
    outLayer = outDataSource.CreateLayer(outputGridfn, spatialRef, geom_type=ogr.wkbPolygon )
    featureDefn = outLayer.GetLayerDefn()

    # create grid cells
    countcols = 0
    while countcols < cols:
        countcols += 1

        # reset envelope for rows
        ringYtop = ringYtopOrigin
        ringYbottom =ringYbottomOrigin
        countrows = 0

        while countrows < rows:
            countrows += 1
            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(ringXleftOrigin, ringYtop)
            ring.AddPoint(ringXrightOrigin, ringYtop)
            ring.AddPoint(ringXrightOrigin, ringYbottom)
            ring.AddPoint(ringXleftOrigin, ringYbottom)
            ring.AddPoint(ringXleftOrigin, ringYtop)
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)

            # add new geom to layer
            outFeature = ogr.Feature(featureDefn)
            outFeature.SetGeometry(poly)
            outLayer.CreateFeature(outFeature)
            outFeature.Destroy

            # new envelope for next poly
            ringYtop = ringYtop - gridHeight
            ringYbottom = ringYbottom - gridHeight

        # new envelope for next poly
        ringXleftOrigin = ringXleftOrigin + gridWidth
        ringXrightOrigin = ringXrightOrigin + gridWidth

    # Close DataSources
    outDataSource.Destroy()


if __name__ == "__main__":

    #
    # example run : $ python grid.py <full-path><output-shapefile-name>.shp xmin xmax ymin ymax gridHeight gridWidth
    #

    #if len( sys.argv ) != 8:
    #    print "[ ERROR ] you must supply seven arguments: output-shapefile-name.shp xmin xmax ymin ymax gridHeight gridWidth"
    #    sys.exit( 1 )
   
    from osgeo import osr

    Proj_CS = osr.SpatialReference()
    Proj_CS.ImportFromEPSG(28355) # This is just an example coordinate system

    xmin='0' 
    xmax='10' 
    ymin='0' 
    ymax='10' 
    gridHeight='1'
    gridWidth='1'
    main(r'C:\Workspace\part0075\MDB modelling\integrated\Modules\Groundwater\GISInterface\GDALInterface\test_scripts_with_GDAL\fish.shp', xmin, xmax, ymin, ymax, gridHeight, gridWidth, Proj_CS)