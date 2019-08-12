# Code derived from:
# https://pcjericks.github.io/py-gdalogr-cookbook/vector_layers.html#create-fishnet-grid

import os
from math import ceil

from osgeo import ogr


def create_fishnet(structured_mesh, spatialRef, copy_dest=None):
    """Function to create a fishnet for the model grid for structured meshes
    Based off of the recipe on:
    https://pcjericks.github.io/py-gdalogr-cookbook/vector_layers.html#create-fishnet-grid

    :param structured_mesh:

    :param spatialRef:

    :param copy_dest: (Default value = None)
    """
    xmin = structured_mesh.xmin
    xmax = structured_mesh.xmax
    ymin = structured_mesh.ymin
    ymax = structured_mesh.ymax
    gridHeight = structured_mesh.gridHeight
    gridWidth = structured_mesh.gridWidth

    outputGridfn = 'structured_model_grid_%im.shp' % int(gridHeight)
    outputGridfnPrj = 'structured_model_grid_%im.prj' % int(gridHeight)

    if copy_dest is not None:
        outputGridDir = os.path.join(copy_dest, 'structured_model_grid_%im' % int(gridHeight))

    # convert sys.argv to float
    xmin = float(xmin)
    xmax = float(xmax)
    ymin = float(ymin)
    ymax = float(ymax)
    gridWidth = float(gridWidth)
    gridHeight = float(gridHeight)

    # get rows
    rows = ceil((ymax - ymin) / gridHeight)
    # get columns
    cols = ceil((xmax - xmin) / gridWidth)

    # start grid cell envelope
    ringXleftOrigin = xmin
    ringXrightOrigin = xmin + gridWidth
    ringYtopOrigin = ymax
    ringYbottomOrigin = ymax - gridHeight

    if not os.path.exists(outputGridDir):
        os.makedirs(outputGridDir)

    # Get the current working directory that was entered into coming into function
    cwd = os.getcwd()
    # Change current working directory to location where grid file will go
    os.chdir(outputGridDir)

    # create output file
    outDriver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(outputGridfn):
        outDriver.DeleteDataSource(outputGridfn)

    outDataSource = outDriver.CreateDataSource(outputGridfn)
    # print dir(outDataSource)
    outLayer = outDataSource.CreateLayer(outputGridfn, spatialRef, geom_type=ogr.wkbPolygon)
    # print dir(outLayer)
    featureDefn = outLayer.GetLayerDefn()
    # print dir(featureDefn)
    # create grid cells
    countcols = 0
    while countcols < cols:
        countcols += 1

        # reset envelope for rows
        ringYtop = ringYtopOrigin
        ringYbottom = ringYbottomOrigin
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

    spatialRef.MorphToESRI()
    with open(outputGridfnPrj, 'w') as f:
        f.write(spatialRef.ExportToWkt())

    # Change back to previous working directory
    os.chdir(cwd)
    outDataSource
    outDataSourceCopy = ogr.GetDriverByName("Memory").CopyDataSource(
        outDataSource, "")

    outDataSource = None  # .Destroy()

    # return outDataSource  # outputGridfn #ds
    return outDataSourceCopy


if __name__ == "__main__":
    # Main to allow standalone testing of this module
    from .GDALInterface import StructuredMesh
    from osgeo import osr

    copy_dest = r"C:\Workspace\part0075\MDB modelling\Campaspe_data\SW\Farm\\"

    Proj_CS = osr.SpatialReference()
    Proj_CS.ImportFromEPSG(28355)  # This is just an example coordinate system

    print(Proj_CS.ExportToWkt())


    structured_mesh = StructuredMesh(xmin='223167.454274',
                                     xmax='223177.454274',
                                     ymin='6052335.91306',
                                     ymax='6052345.91306',
                                     gridHeight='1',
                                     gridWidth='1')

    test = create_fishnet(structured_mesh, Proj_CS, copy_dest=copy_dest)
    #test = None
