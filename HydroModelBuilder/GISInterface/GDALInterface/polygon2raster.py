import numpy as np
from osgeo import gdal, ogr

RASTERIZE_COLOR_FIELD = "__color__"


def test_raster(raster_obj):
    """Function to test if raster generated has all of the features
    defined in the polygon from which it was generated.

    :param raster_obj:
    """
    raster_array = raster_obj.ReadAsArray()
    return max([len(np.unique(raster_array[i])) - 1 for i in range(3)])
# End test_raster()


def rasterize(source_ds, out_fname='default.tif', pixel_size=None, bounds=None, feature_name=None,
              field_type=ogr.OFTInteger):
    """Function to create a raster based on a polygon shapefile. If it is a multi-
    polygon shapefile, then a feature name can be used to split the raster up.


    :param source_ds:
        :param out_fname:  (Default value = 'default.tif')
    :param pixel_size: (Default value = None)
    :param bounds: (Default value = None)
    :param feature_name: (Default value = None)
    :param field_type: (Default value = ogr.OFTInteger)
    :param out_fname:  (Default value = 'default.tif')
    """
    feature_dict = {}

    if pixel_size == None:
        pixel_size = 25

    source_layer = source_ds.GetLayer(0)

    source_srs = source_layer.GetSpatialRef()
    if bounds == None:
        x_min, x_max, y_min, y_max = source_layer.GetExtent()
    else:
        x_min, x_max, y_min, y_max = bounds

    # Create a field in the source layer to hold the features colors
    field_def = ogr.FieldDefn(RASTERIZE_COLOR_FIELD, field_type)
    source_layer.CreateField(field_def)
    source_layer_def = source_layer.GetLayerDefn()
    field_index = source_layer_def.GetFieldIndex(RASTERIZE_COLOR_FIELD)
    # Generate random values for the color field (it's here that the value
    # of the attribute should be used, but you get the idea)
    for index, feature in enumerate(source_layer):
        idx = index + 1
        feature.SetField(field_index, idx)
        source_layer.SetFeature(feature)
        if feature_name is not None:
            feature_dict[idx] = feature.GetField(feature_name)
    # Create the destination data source
    x_res = int((x_max - x_min) / pixel_size)
    y_res = int((y_max - y_min) / pixel_size)

    target_ds = gdal.GetDriverByName('GTiff').Create(out_fname + '.tif', x_res,
                                                     y_res, 3, gdal.GDT_Byte)
    target_ds.SetGeoTransform((
        x_min, pixel_size, 0,
        y_max, 0, -pixel_size,
    ))

    if source_srs:
        # Make the target raster have the same projection as the source
        target_ds.SetProjection(source_srs.ExportToWkt())
    else:
        # Source has no projection (needs GDAL >= 1.7.0 to work)
        target_ds.SetProjection('LOCAL_CS["arbitrary"]')
    # Rasterize
    err = gdal.RasterizeLayer(target_ds, (3, 2, 1), source_layer,
                              burn_values=(0, 0, 0),
                              options=["ATTRIBUTE=%s" % RASTERIZE_COLOR_FIELD, "ALL_TOUCHED=FALSE"])
    if err != 0:
        raise Exception("error rasterizing layer: %s" % err)

    if source_layer.GetFeatureCount() > test_raster(target_ds):
        print("Warning: Unique values in raster less than feature count !!!")
    # End if

    # return target_ds_copy, feature_dict
    return target_ds, feature_dict
# End rasterize()


def array_from_rasterize(source_ds, out_fname='default.tif', pixel_size=None,
                         bounds=None, feature_name=None, field_type=ogr.OFTInteger):
    """

    :param source_ds:
    :param out_fname:  (Default value = 'default.tif')
    :param pixel_size: (Default value = None)
    :param bounds: (Default value = None)
    :param feature_name: (Default value = None)
    :param field_type: (Default value = ogr.OFTInteger)
    """
    raster, feature_dict = rasterize(source_ds, out_fname=out_fname,
                                     pixel_size=pixel_size, bounds=bounds,
                                     feature_name=feature_name, field_type=ogr.OFTInteger)

    return raster.ReadAsArray()[0], feature_dict
