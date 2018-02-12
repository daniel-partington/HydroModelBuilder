# Taken from http://stackoverflow.com/questions/15968762/shapefile-and-matplotlib-plot-polygon-collection-of-shapefile-coordinates

import matplotlib.patches as mpatches
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import numpy as np
from osgeo import ogr


def plotShpfile(shpfile=None, shptype=None, facecolor='gray', edgecolor='black'):
    """

    :param shpfile: (Default value = None)
    :param shptype: (Default value = None)
    :param facecolor: (Default value = 'gray')
    :param edgecolor: (Default value = 'black')
    """
    if type(shpfile) == str:
        shpfile = [shpfile]
        shptype = [shptype]

    for index, shp in enumerate(shpfile):
        # Extract first layer of features from shapefile using OGR
        ds = None

        ds = ogr.Open(shpfile[index])
        nlay = ds.GetLayerCount()
        lyr = ds.GetLayer(0)

        if index == 0:
            # Get extent and calculate buffer size
            ext = lyr.GetExtent()
            xoff = (ext[1] - ext[0]) / 50
            yoff = (ext[3] - ext[2]) / 50

            # Prepare figure
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_xlim(ext[0] - xoff, ext[1] + xoff)
            ax.set_ylim(ext[2] - yoff, ext[3] + yoff)

        shptype2 = shptype[index]
        if shptype2 == 'points':
            for feat in lyr:
                geomRef = feat.GetGeometryRef()
                x = [geomRef.GetX(i) for i in range(geomRef.GetPointCount())]  # although points has only one object..
                y = [geomRef.GetY(i) for i in range(geomRef.GetPointCount())]  # although points has only one object..
                i1, = ax.plot(x, y, 'ro', alpha=0.5)
            # ax.set_aspect(1.0)
            # plt.show()

        if shptype2 == 'polyline':
            for feat in lyr:
                geomRef = feat.GetGeometryRef()
                x = [geomRef.GetX(i) for i in range(geomRef.GetPointCount())]
                y = [geomRef.GetY(i) for i in range(geomRef.GetPointCount())]
                i2, = ax.plot(x, y, color='blue')
            # ax.set_aspect(1.0)
            # plt.show()

        if shptype2 == 'polygon':
            paths = []
            lyr.ResetReading()

            # Read all features in layer and store as paths
            for feat in lyr:
                geom = feat.geometry()
                codes = []
                all_x = []
                all_y = []
                for i in range(geom.GetGeometryCount()):
                    # Read ring geometry and create path
                    r = geom.GetGeometryRef(i)
                    x = [r.GetX(j) for j in range(r.GetPointCount())]
                    y = [r.GetY(j) for j in range(r.GetPointCount())]
                    # skip boundary between individual rings
                    codes += [mpath.Path.MOVETO] + \
                        (len(x) - 1) * [mpath.Path.LINETO]
                    all_x += x
                    all_y += y
                path = mpath.Path(np.column_stack((all_x, all_y)), codes)
                paths.append(path)

            # Add paths as patches to axes
            for path in paths:
                patch = mpatches.PathPatch(path,
                                           facecolor=facecolor, edgecolor=edgecolor, alpha=0.5)
                ax.add_patch(patch)

    ax.set_aspect(1.0)
    plt.show()


if __name__ == "__main__":
    shpfile0 = r"C:\Workspace\part0075\MDB modelling\testbox\01_steady_state\GW_model_area_model.shp"
    #plotShpfile(shpfile=shpfile, shptype='polygon')
    shpfile1 = r"C:\Workspace\part0075\MDB modelling\testbox\01_steady_state\Campaspe_Riv_model.shp"
    #plotShpfile(shpfile=shpfile, shptype='polyline')
    shpfile2 = r"C:\Workspace\part0075\MDB modelling\testbox\01_steady_state\pumping wells_reproj.shp"

    plotShpfile(shpfile=[shpfile0, shpfile1, shpfile2], shptype=['polygon', 'polyline', 'points'])
