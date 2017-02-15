# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 08:11:11 2016

@author: part0075

Mesh builder from stratigraphy

Build from top to bottom.

This is only to be called once the rasters have been reprojected into the
correct coordinate system

"""
import os

import matplotlib.pyplot as plt
import numpy as np
from osgeo import gdal, gdalconst
from scipy import ndimage as nd

import array2Vtk


# Step 1. Load raster layers top and bottom


def map_raster_array_to_mesh(hu_raster_path, hu_raster_files, out_path, vtk_out,
                             min_height, max_height):
    raster_set = {}

    print hu_raster_path

    for raster in hu_raster_files:
        fname = os.path.join(hu_raster_path, raster) # + hu_ext
        print 'Processing: ', fname    
        ds = gdal.Open(fname, gdalconst.GA_ReadOnly)    
        raster_set[raster] = [ds.GetRasterBand(1).ReadAsArray(), 0.0] 

        ds = gdal.Open(fname, gdalconst.GA_ReadOnly)
        raster_set[raster] = [ds.GetRasterBand(1).ReadAsArray(), 0.0]

        if raster == hu_raster_files[0]:
            ncol = ds.RasterXSize
            nrow = ds.RasterYSize
            (x0, delr, z0, y0, z1, delc) = ds.GetGeoTransform()
            delc = -delc
        ds = None

    # This solution taken from:
    # http://stackoverflow.com/questions/3662361/fill-in-missing-values-with-nearest-neighbour-in-python-numpy-masked-arrays
    def fill(data, invalid=None):
        """
        Replace the value of invalid 'data' cells (indicated by 'invalid')
        by the value of the nearest valid data cell

        Input:
            data:    numpy array of any dimension
            invalid: a logical array of same shape as 'data'. True cells set 
                     where data value should be replaced.
                     If None (default), use: invalid  = np.isnan(data)

        Output:
            Return a filled array.
        """

        if invalid is None:
            invalid = np.isnan(data)

        ind = nd.distance_transform_edt(invalid, return_distances=False, return_indices=True)
        return data[tuple(ind)]

    for raster in hu_raster_files:
        data = raster_set[raster][0]
        invalid = np.ma.masked_array(data, data == raster_set[raster][1])
        raster_set[raster] += [invalid]
        a = fill(data, invalid=invalid.mask)
        raster_set[raster] += [a]

        # Clean_up
        data = None
        invalid = None
        a = None
    # end for

    # Determine required number of layers based on
    # Note: this is integer division
    layers = len(hu_raster_files) / 2

    mesh = np.zeros((layers + 1,
                     len(raster_set[hu_raster_files[0]][0]),
                     len(raster_set[hu_raster_files[0]][0][0])))

    layer_set = (layers, len(raster_set[hu_raster_files[0]][0]),
                 len(raster_set[hu_raster_files[0]][0][0]))

    zone_matrix = np.zeros(layer_set)
    thickness = np.zeros(layer_set)

    # Create an array of ignored and applied values to keep track of processing
    mask_active = np.full(layer_set, False, dtype=bool)

    mask_thickness = np.full(layer_set, False, dtype=bool)

    cumulative_active = np.full(layer_set, False, dtype=bool)

    for index, raster in enumerate(hu_raster_files):
        if index % 2 == 1:
            continue
        mask_active[index / 2] = raster_set[raster][2].mask

        # Test if any of the cells are too thin and set to inactive if so, but
        # maintain the top of the layer to prevent mismatches between layer top
        # and bottoms

        mask_thickness[index / 2] = raster_set[raster][3] - \
            raster_set[hu_raster_files[index + 1]][3] < min_height
        mask_active[index / 2] = mask_active[index / 2] | mask_thickness[index / 2]
    # End for

    index = 0
    # Set top
    mesh[index / 2] = raster_set[hu_raster_files[index]][3]
    # Set bottom of top layer
    mesh[index / 2 + 1] = raster_set[hu_raster_files[index + 1]][3]

    zone_matrix[index / 2][mask_active[index / 2]] = -1
    zone_matrix[index / 2][~mask_active[index / 2]] = index / 2 + 1

    for i in xrange(len(hu_raster_files)):
        if i % 2 == 1:
            continue
        cumulative_active[i / 2] = ~mask_active[i / 2]
    # end for

    for i in xrange(index + 2, len(hu_raster_files)):
        if i % 2 == 1:
            continue

        # Modify properties of layer using merge up logical array
        for j in xrange(0, i):
            if j % 2 == 1:
                continue

            # Identify intersection of where this layer is inactive and where
            # layer below is active
            merge_up_locations = ~(cumulative_active[j / 2] | mask_active[i / 2])
            # Modify top of mesh in this layer to the merge up locations,
            # but excluding those where the thickness mask was applied

            # Find where existing mesh is higher than candidate replacement from
            # next layer, so we know where to keep existing mesh
            new_elev = raster_set[hu_raster_files[i]][3] - mesh[j / 2] > 0

            # Find cells where they need updating and where cells aren't already
            # defined by zone, i.e. inactive cells
            modify_cells = new_elev & ~(cumulative_active[j / 2] | ~mask_active[i / 2])

            if j > 0:
                # Check above layer which might have active cells for which we don't want
                # to modify the bottom of the layer
                combined = (merge_up_locations | modify_cells) & ~cumulative_active[j / 2 - 1]
            else:
                # Combine the two above logical arrays to apply in mesh modification
                combined = merge_up_locations | modify_cells

            # ----------------------------------------------------------------------
            mesh[j / 2][combined] = raster_set[hu_raster_files[i]][3][combined]
            # ----------------------------------------------------------------------

            if j > 0:
                test = mesh[j / 2 - 1] - mesh[j / 2] < 0
                # if np.any(test == True):
                #    print 'Error when processing layer test2', i/2

            zone_matrix[j / 2][merge_up_locations] = i / 2 + 1

            # Modify the intermediate mesh levels for those layers that are
            # split over multiple layers to retain original thickness

            # print 'Working on layer: ', i/2
            # print 'Modifying mesh: ', j/2+1
            # print 'Multiplier = ', (j/2+1)
            # print 'Divisor = ', (i/2+1)

            proposed_heights = (raster_set[hu_raster_files[i]][3][merge_up_locations] -
                                (j / 2 + 1) * ((raster_set[hu_raster_files[i]][3]
                                                [merge_up_locations] -
                                                raster_set[hu_raster_files[i + 1]
                                                           ][3][merge_up_locations]) /
                                               (i / 2 + 1)))

            mesh[j / 2 + 1][merge_up_locations] = proposed_heights

            test = mesh[j / 2] - mesh[j / 2 + 1] < 0
            # if np.any(test == True):
            #    print 'Error when processing layer test3', i/2

            # After modifying the layer, update cumulative active for modified layer
            cumulative_active[j / 2] = cumulative_active[j / 2] | ~mask_active[i / 2]

        # Set properties for layer from i/2
        zone_matrix[i / 2][mask_active[i / 2]] = -1
        zone_matrix[i / 2][~mask_active[i / 2]] = i / 2 + 1
        # Set the base for layer i/2

        # Implement bottom of mesh for this layer BUT first check that the
        # inactive cell levels are not higher than the tpo of layer.
        # First identify where such cells might occur
        level_test = raster_set[hu_raster_files[i + 1]][3] - mesh[i / 2] < 0
        # Get array of inactive cells in layer:
        # layer_inactive = ~cumulative_active[i/2]
        # Identify logical array of cells to modify:
        # modify = level_test & layer_inactive

        mesh[i / 2 + 1][level_test] = raster_set[hu_raster_files[i + 1]][3][level_test]

    # End for

    # Use the fill function to clean up the inactive parts of the mesh to assist
    # with the rendering

    for index, raster in enumerate(hu_raster_files):
        if index % 2 == 1:
            continue
        if index / 2 == 0:
            mesh[index / 2] = fill(mesh[index / 2], invalid=np.ma.masked_array(
                zone_matrix[index / 2], zone_matrix[index / 2] == -1).mask)
        elif index < len(hu_raster_files):
            mesh[index / 2] = fill(mesh[index / 2],
                                   invalid=np.ma.masked_array(zone_matrix[index / 2],
                                                              zone_matrix[index / 2] == -1).mask |
                                   np.ma.masked_array(zone_matrix[index / 2 - 1],
                                                      zone_matrix[index / 2 - 1] == -1).mask)

        if index / 2 == len(hu_raster_files) / 2 - 1:
            mesh[index / 2 + 1] = fill(mesh[index / 2 + 1], invalid=np.ma.masked_array(
                zone_matrix[index / 2], zone_matrix[index / 2] == -1).mask)

    # Check that all layers are not intersecting
    for index, raster in enumerate(hu_raster_files):
        if index % 2 == 1:
            continue
        test = mesh[index / 2] - mesh[index / 2 + 1] <= 1
        if np.any(test == True):
            print 'Error when processing layer ', index / 2

    # Calculate the thickness of each layer on the final mesh
    for index, raster in enumerate(hu_raster_files):
        # Only work on even index
        if index % 2 == 1:
            continue

        thickness[index / 2] = mesh[index / 2] - mesh[index / 2 + 1]
    # End for

    tester2 = True
    if tester2 == True:
        raster_thickness = {}

        for index, raster in enumerate(hu_raster_files):
            if index % 2 == 1:
                continue
            raster_thickness[index / 2] = np.ma.masked_where(
                (zone_matrix[index / 2] != index / 2 + 1), raster_set[raster][3] -
                raster_set[hu_raster_files[index + 1]][3])
            print 'Using', raster, ' and ', hu_raster_files[index + 1]

        mesh_zone_thickness = {}

        for i in xrange(len(hu_raster_files) / 2):
            top = np.zeros((len(raster_set[hu_raster_files[0]][0]),
                            len(raster_set[hu_raster_files[0]][0][0])))
            bot = np.zeros((len(raster_set[hu_raster_files[0]][0]),
                            len(raster_set[hu_raster_files[0]][0][0])))

            bot[zone_matrix[i] == i + 1] = mesh[i + 1][zone_matrix[i] == i + 1]

            for layer in reversed(range(zone_matrix.shape[0])):
                top[zone_matrix[layer] == i + 1] = mesh[layer][zone_matrix[layer] == i + 1]

            mesh_zone_thickness[i] = np.ma.masked_where((zone_matrix[i] != i + 1), top - bot)

        for i in xrange(thickness.shape[0]):
            fig = plt.figure()
            fig.add_subplot(1, 3, 1, aspect='equal')        
            plt.imshow(mesh_zone_thickness[i], interpolation='none')
            plt.title('Thickness in mesh: ' + hu_raster_files[i * 2])
            plt.colorbar()
            fig.add_subplot(1, 3, 2, aspect='equal')        
            plt.imshow(raster_thickness[i], interpolation='none')        
            plt.title('Thickness in raster')
            plt.colorbar()
            fig.add_subplot(1, 3, 3, aspect='equal')        
            plt.imshow(mesh_zone_thickness[i] - raster_thickness[i], interpolation='none')        
            plt.title('Difference in thickness')
            plt.colorbar()

    tester = False
    if tester == True:
        # print thickness
        thickness2 = {}
        thick_zero = {}

        thick_shape = thickness.shape[0]

        for i in xrange(thick_shape):
            thickness2[i] = np.ma.masked_where((zone_matrix[i] < 0), thickness[i])
            thick_zero[i] = np.ma.masked_where((zone_matrix[i] != 0), thickness[i])

        for i in xrange(thick_shape):
            fig = plt.figure()
            fig.add_subplot(1, 3, 1, aspect='equal')        
            plt.imshow(thickness2[i], interpolation='none')
            plt.title('Thickness in active cells: ', hu_raster_files[i * 2][0:5])
            plt.colorbar()
            fig.add_subplot(1, 3, 2, aspect='equal')        
            plt.imshow(thick_zero[i], interpolation='none')        
            plt.title('Areas where thickness is 0')
            plt.colorbar()
            fig.add_subplot(1, 3, 3, aspect='equal')        
            plt.imshow(zone_matrix[i], interpolation='none')        
            plt.title('Zonal delineation')
            plt.colorbar()

    # for i in range(zone_matrix.shape[0]):
    #    plt.figure()
    #    plt.imshow(zone_matrix[i], interpolation='none')

    if np.any(thickness < 0.):
        print 'issues'

    # for i in range(mesh.shape[0]):
    #    plt.figure()
    #    plt.imshow(mesh[i], interpolation='none')

    grid_info = [ncol, nrow, delc, delr, x0, y0]
    array2Vtk.build_vtk_from_array(grid_info, np.fliplr(mesh), ["z_elev"], [np.fliplr(mesh)], [
                                   "zone", "thickness"],
                                   [np.fliplr(zone_matrix), np.fliplr(thickness)],
                                   out_path, vtk_out)

    return mesh, zone_matrix, grid_info


if __name__ == "__main__":
    bottom_ext = "_2b"
    top_ext = "_1t"

    # hu_raster_path = r"C:\Workspace\part0075\MDB modelling\VAF_v2.0_ESRI_GRID\ESRI_GRID\Preprocessed_data\Hydrogeological_Unit_Layers\\"
    # hu_raster_files = ["qa_1t_bb", "qa_2b_bb", "utqa_1t_bb", "utqa_2b_bb",
    # "utaf_1t_bb", "utaf_2b_bb", "lta_1t_bb", "lta_2b_bb", "cps_1t_bb",
    # "cps_2b_bb"]

    hu_raster_path = r"C:\Workspace\part0075\MDB modelling\integrated\Modules\Groundwater\model_files\structured_model_grid_\\"
    hu_raster_files = ["qa_1t_model_grid.bil", "qa_2b_model_grid.bil", "utb_1t_model_grid.bil", "utb_2b_model_grid.bil", "utqa_1t_model_grid.bil", "utqa_2b_model_grid.bil", "utam_1t_model_grid.bil",
                       "utam_2b_model_grid.bil", "utaf_1t_model_grid.bil", "utaf_2b_model_grid.bil", "lta_1t_model_grid.bil", "lta_2b_model_grid.bil", "bse_1t_model_grid.bil", "bse_2b.tif_model_grid.bil"]

    dx = "1000m"

    hu_raster_path = hu_raster_path[:-2] + dx + hu_raster_path[-2:]
    # hu_raster_files = [x[:-4] + '_' + dx + x[-4:] for x in hu_raster_files]

    #      "utb_1t_model_grid.bil", "utb_2b_model_grid.bil",
    #      "qa_1t_model_grid.bil", "qa_2b_model_grid.bil",   , "bse_1t_model_grid.bil", "bse_2b.tif_model_grid.bil"

    # hu_layers = {"qa":["qa_1t_mb", "qa_2b_mb"], "utb":["utb_1t_mb", "utb_2b_mb"], "utqa":["utqa_1t_mb", "utqa_2b_mb"], "utaf":["utaf_1t_mb", "utaf_2b_mb"], "lta":["lta_1t_mb", "lta_2b_mb"], "bse":["bse_1t_mb", "bse_2b_mb.tif"]}
    hu_layers_ordering = ["qa", "utb", "utqa", "utaf", "lta", "bse"]
    # hu_ext = "_mb"
    vtk_out = 'Campaspe_all_model_mesh'

    out_path = r"C:\Workspace\part0075\MDB modelling\integrated\Modules\Groundwater\model_files\testing\\"

    # Minimum thickness of cells
    min_height = 1
    # Maximum thickness of cells
    max_height = 1000

    map_raster_array_to_mesh(hu_raster_path, hu_raster_files,
                             out_path, vtk_out, min_height, max_height)
