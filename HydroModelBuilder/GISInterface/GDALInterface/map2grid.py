import os
import sys

from osgeo import ogr


def shp2grid(shp_to_map, poly_mesh, shp_type=None, feature_id=None, data_folder=None):
    """

    :param shp_to_map:
    :param poly_mesh:
    :param shp_type: (Default value = None)
    :param feature_id: (Default value = None)
    :param data_folder: (Default value = None)
    """
    pwd = os.getcwd()
    if data_folder != None:
        os.chdir(data_folder)

    if shp_type in ['poly', 'points']:
        print 'Processing shapefile with ', shp_type
    else:
        print 'Shape type not recognised: ', shp_type
        return
    # end if

    mesh_layer = poly_mesh.GetLayer()
    srs = mesh_layer.GetSpatialRef()

    shp_to_map_layer = shp_to_map.GetLayer()

    # Set driver for creating shape files
    driver_inter = ogr.GetDriverByName('MEMORY')

    shp_to_map_mapped = []

    mesh_layer.ResetReading()

    x = mesh_layer.GetFeatureCount()

    for feature in mesh_layer:

        # Track progress because this is very slow
        i = feature.GetFID()
        sys.stdout.write('\r')
        sys.stdout.write("[%-20s] %d%%" %
                         ('=' * int(float(i) / float(x) * 20. + 1), int(100.0 / (x - 1.0) * i)))
        sys.stdout.flush()

        grid_cell = feature.GetGeometryRef()

        mesh_layer_ds = None
        while mesh_layer_ds == None:
            #mesh_layer_ds = driver_inter.CreateDataSource('mesh_temp.shp')
            mesh_layer_ds = driver_inter.CreateDataSource('')

        mesh_temp = None
        while mesh_temp == None:
            mesh_temp = mesh_layer_ds.CreateLayer('mylayer', srs, geom_type=ogr.wkbMultiPolygon)

        featureDefn = mesh_temp.GetLayerDefn()

        ring = ogr.Geometry(ogr.wkbLinearRing)
        (xmin, xmax, ymin, ymax) = grid_cell.GetEnvelope()
        ring.AddPoint(xmin, ymax)
        ring.AddPoint(xmax, ymax)
        ring.AddPoint(xmax, ymin)
        ring.AddPoint(xmin, ymin)
        ring.AddPoint(xmin, ymax)

        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)
        # add new geom to layer
        outFeature = ogr.Feature(featureDefn)
        outFeature.SetGeometry(poly)
        mesh_temp.CreateFeature(outFeature)
        outFeature.Destroy

        dstshp = None
        while dstshp == None:
            #dstshp = driver_inter.CreateDataSource('temp.shp')
            dstshp = driver_inter.CreateDataSource('')

        if shp_type == 'points':
            dstlayer = dstshp.CreateLayer('mylayer', srs, geom_type=ogr.wkbPoint)
        elif shp_type == 'poly':
            dstlayer = dstshp.CreateLayer('mylayer', srs, geom_type=ogr.wkbLineString)

        shp_to_map_layer.Intersection(mesh_temp, dstlayer)

        if int(dstlayer.GetFeatureCount()) != 0:
            if shp_type == 'points':
                point_ids = []
                for feature2 in dstlayer:
                    point_ids += [feature2.GetField(feature_id)]
                # end for
                centroid_txt = grid_cell.Centroid().ExportToWkt()
                centroid_nums = centroid_txt.split('(')[1].split(')')[0]
                centroid = centroid_nums.split(' ')
                shp_to_map_mapped += [[point_ids, centroid]]
            elif shp_type == 'poly':
                length = 0.0
                for feature2 in dstlayer:
                    line = feature2.GetGeometryRef()
                    length += line.Length()
                # print length, 'm'
                # print grid_cell.Centroid()
                centroid_txt = grid_cell.Centroid().ExportToWkt()
                centroid_nums = centroid_txt.split('(')[1].split(')')[0]
                centroid = centroid_nums.split(' ')
                shp_to_map_mapped += [[length, centroid]]
            # end if
        # end if

        # clean up
        dstlayer = None
        dstshp = None
        mesh_temp = None
        mesh_layer_ds = None

    # close shape files
    mesh_layer = None
    shp_to_map_layer = None

    # Go back to working directory
    os.chdir(pwd)

    return shp_to_map_mapped


if __name__ == "__main__":
    # Open a points object
    #    driver = ogr.GetDriverByName("ESRI Shapefile")
    #    ds = driver.Open(r"C:\Workspace\part0075\MDB modelling\testbox\data_build\pumping wells_reproj.shp", 0)
    #    poly_obj = ds.GetLayer()
    #    if poly_obj == None:
    #        print 'Could not open '
    #    srs = poly_obj.GetSpatialRef()
    #
    #    # Open the mesh object
    #    #ds2 = driver.Open(r"C:\Workspace\part0075\MDB modelling\integrated\Modules\Groundwater\model_files\structured_model_grid_1000m\structured_model_grid_1000m.shp", 0)
    #    ds2 = driver.Open(r"C:\Workspace\part0075\MDB modelling\Campaspe_data\SW\Farm\structured_model_grid_5000m.shp", 0)
    #
    #    mapped_list = shp2grid(ds, ds2, shp_type='points', feature_id="OLD ID")
    #    # for item in mapped_list:
    #    #    print item
    #    print mapped_list  # [0]
    #
    #    ds = None
    #    ds2 = None

    # Open a polyline object
    driver = ogr.GetDriverByName("ESRI Shapefile")
    ds = driver.Open(r"C:\Workspace\part0075\MDB modelling\Campaspe_data\SW\Farm\test_model.shp", 0)
    #ds = driver.Open(r"C:\Workspace\part0075\MDB modelling\testbox\input_data\Waterways\Campaspe_Riv.shp", 0)
    poly_obj = ds.GetLayer()
    if poly_obj == None:
        print 'Could not open '
    srs = poly_obj.GetSpatialRef()

    # Open the mesh object
    ds2 = driver.Open(
        r"C:\Workspace\part0075\MDB modelling\testbox\00_Campaspe_Cascade\02_transient_flow\structured_model_grid_1000m\structured_model_grid_1000m.shp", 0)
    #ds2 = driver.Open(r"C:\Workspace\part0075\MDB modelling\Campaspe_data\SW\Farm\structured_model_grid_5000m.shp", 0)

    import time

    start1 = time.time()
    mapped_list = shp2grid(ds, ds2, shp_type='poly')
    end1 = time.time()
    # for item in mapped_list:
    #    print item

    ds = None
    ds2 = None

    # Open a polyline object
    driver = ogr.GetDriverByName("ESRI Shapefile")
    ds = driver.Open(r"C:\Workspace\part0075\MDB modelling\Campaspe_data\SW\Farm\test_model_new.shp", 0)
    poly_obj = ds.GetLayer()
    if poly_obj == None:
        print 'Could not open '
    srs = poly_obj.GetSpatialRef()

    # Open the mesh object
    ds2 = driver.Open(
        r"C:\Workspace\part0075\MDB modelling\testbox\00_Campaspe_Cascade\02_transient_flow\structured_model_grid_1000m\structured_model_grid_1000m.shp", 0)
    #ds2 = driver.Open(r"C:\Workspace\part0075\MDB modelling\Campaspe_data\SW\Farm\structured_model_grid_5000m.shp", 0)

    start2 = time.time()
    mapped_list = shp2grid(ds, ds2, shp_type='poly')
    end2 = time.time()
    # for item in mapped_list:
    #    print item

    print(end1 - start1, end2 - start2)

    ds = None
    ds2 = None
