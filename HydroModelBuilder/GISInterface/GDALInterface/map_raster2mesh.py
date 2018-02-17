# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 08:11:11 2016

@author: part0075

Mesh builder from stratigraphy

Build from top to bottom.

This is only to be called once the rasters have been reprojected into the
correct coordinate system

"""
import itertools as it
import os
from itertools import groupby

import matplotlib.pyplot as plt
import numpy as np
from osgeo import gdal, gdalconst
from scipy import ndimage as nd

import array2Vtk


def reclassIsolatedCells(mesh3D_1, passes=1, assimilate=False):
    """Remove cells that are surrounded by non-active cells in the horizontal plane

    e.g. if cells with positive integer is surrounded in above, below and to each side, then reassign to -1.

    :param mesh3D_1: ndarray, array to remove isolated cells from
    :param passes: int, number of times to go over the array (Default value = 1)
    :param assimilate: ? Default value = False)
    """
    # Clean up idle cells:
    (lay, row, col) = mesh3D_1.shape
    loop_elements = it.product(range(passes), range(lay), range(row), range(col))
    for p, k, j, i in loop_elements:
        cell_zone = mesh3D_1[k][j][i]
        if cell_zone == -1:
            continue
        # End if

        target_zone = mesh3D_1[k]
        # Assimilate cell if surrounded by four of the same
        if assimilate:
            if (((j > 0) and (target_zone[j - 1][i] != cell_zone)) and        # North
                ((j < row - 1) and (target_zone[j + 1][i] != cell_zone)) and  # South
                ((i < col - 1) and (target_zone[j][i + 1] != cell_zone)) and  # East
                ((i > 0) and (target_zone[j][i - 1] != cell_zone)) and        # West
                # North-East
                ((j > 0) and (i < col - 1) and (target_zone[j - 1][i + 1] != cell_zone)) and
                # South-East
                ((j < row - 1) and (i < col - 1) and (target_zone[j + 1][i + 1] != cell_zone)) and
                # North-West
                ((j > 0) and (i > 0) and (target_zone[j - 1][i - 1] != cell_zone)) and
                    ((j < row - 1) and (i > 0) and (target_zone[j + 1][i - 1] != cell_zone))):  # South-West

                neighbours = []
                if j > 0:
                    neighbours += [target_zone[j - 1][i]]
                if j < row - 1:
                    neighbours += [target_zone[j + 1][i]]
                if i < col - 1:
                    neighbours += [target_zone[j][i + 1]]
                if i > 0:
                    neighbours += [target_zone[j][i - 1]]
                # End if

                def most_common_oneliner(L):
                    """
                    :param L:
                    """
                    return max(groupby(sorted(L)), key=lambda(x, v): (len(list(v)), -L.index(x)))[0]

                most_common = most_common_oneliner(neighbours)
                if most_common != -1:
                    target_zone[j][i] = most_common_oneliner(neighbours)

        # Check North, South, East, West zones
        # If any condition is true, then continue on
        # Otherwise, set the cell to -1
        if (((j > 0) and (target_zone[j - 1][i] != -1)) or        # North
            ((j < row - 1) and (target_zone[j + 1][i] != -1)) or  # South
            ((i < col - 1) and (target_zone[j][i + 1] != -1)) or  # East
                ((i > 0) and (target_zone[j][i - 1] != -1))):     # West
            continue
        # End if

        # None of the above conditions were true
        target_zone[j][i] = -1
        print("reclassifying")
    # End for
    return mesh3D_1
# End reclassIsolatedCells()


def fill(data, invalid=None):
    """Replace the value of invalid 'data' cells with the value of the nearest valid data cell.

    This solution taken from:
    http://stackoverflow.com/questions/3662361/fill-in-missing-values-with-nearest-neighbour-in-python-numpy-masked-arrays

    :param data: ndarray, raster
    :param invalid: ndarray, boolean array of same shape as `data`.
                    `True` indicates cells to be replaced.
                    If `None` (default), invalid cells will be determined with `np.isnan(data)`

    :returns: ndarray
    """
    if invalid is None:
        invalid = np.isnan(data)

    ind = nd.distance_transform_edt(invalid, return_distances=False, return_indices=True)
    return data[tuple(ind)]
# End fill()


def map_raster_array_to_mesh(hu_raster_path, hu_raster_files, out_path, vtk_out,
                             min_height, max_height):
    """

    :param hu_raster_path: str, path to files
    :param hu_raster_files: list[str], hydrological unit raster files
    :param out_path: str, output path
    :param vtk_out: str, name to use for VTK output
    :param min_height: float, minimum cell height
    :param max_height: float, maximum cell height (not used?)
    """
    raster_set = {}

    print("HU Raster Path:", hu_raster_path)
    for raster in hu_raster_files:
        fname = os.path.join(hu_raster_path, raster)
        # print('Processing: ', fname)
        ds = gdal.Open(fname, gdalconst.GA_ReadOnly)

        raster_set[raster] = {
            "original": ds.GetRasterBand(1).ReadAsArray(),
            "no_data": ds.GetRasterBand(1).GetNoDataValue()
        }

        if raster == hu_raster_files[0]:
            # Get raster info from first file as these should be identical
            ncol = ds.RasterXSize
            nrow = ds.RasterYSize
            (x0, delr, z0, y0, z1, delc) = ds.GetGeoTransform()
            delc = -delc
        # End if

        ds = None  # destroy GDAL object
    # End for

    for raster in hu_raster_files:
        data = raster_set[raster]['original']
        invalid = np.ma.masked_equal(data, np.min(data))
        raster_set[raster]['invalid'] = invalid
        raster_set[raster]['filled'] = fill(data, invalid=invalid.mask)
    # End for

    # Clean_up
    data = None
    invalid = None

    # Determine required number of layers based on
    # Note: this is integer division
    num_hu_files = len(hu_raster_files)
    num_layers = num_hu_files / 2

    # Map raster id, layer id, and name of raster file [i, int(i / 2), filename]
    raster_id_map = zip([i for i in range(num_hu_files + 1)], [i / 2 for i in range(num_hu_files + 1)], hu_raster_files)

    top_most_raster = raster_set[hu_raster_files[0]]['original']
    mesh = np.zeros((num_layers + 1,
                     len(top_most_raster),
                     len(top_most_raster[0])))

    layer_set = (num_layers, len(top_most_raster), len(top_most_raster[0]))
    zone_matrix = np.zeros(layer_set)
    thickness = np.zeros(layer_set)

    # Create an array of ignored and applied values to keep track of processing
    inactive_cells = np.full(layer_set, False, dtype=bool)  # keeps track of inactive cells
    mask_thickness = np.full(layer_set, False, dtype=bool)  # keep track of layer thickness
    cumulative_active = np.full(layer_set, False, dtype=bool)  # keeps track of active cells

    top_rasters = raster_id_map[::2]
    for idx, mod_idx, raster in top_rasters:
        this_raster = raster_set[raster]['filled']
        next_raster = raster_set[hu_raster_files[idx + 1]]['filled']
        inactive_cells[mod_idx] = raster_set[raster]['invalid'].mask

        # Test if any of the cells are too thin and set to inactive if so, but
        # maintain the top of the layer to prevent mismatches between layer top
        # and bottoms
        mask_thickness[mod_idx] = (this_raster - next_raster) < min_height
        inactive_cells[mod_idx] = inactive_cells[mod_idx] | mask_thickness[mod_idx]
    # End for

    mesh[0] = raster_set[hu_raster_files[0]]['filled']  # Set top
    mesh[1] = raster_set[hu_raster_files[1]]['filled']  # Set bottom of top layer

    zone_matrix[0][inactive_cells[0]] = -1  # Set inactive cells for zone to -1
    zone_matrix[0][~inactive_cells[0]] = 1  # Set Zone IDs

    for i in xrange(1, num_hu_files, 2):
        mod_idx = i / 2
        cumulative_active[mod_idx] = ~inactive_cells[mod_idx]
    # End for

    # loop over 'top' rasters, skipping rasters for the first layer (rasters 0 and 1)
    for i in xrange(2, num_hu_files, 2):
        mod_i_idx = i / 2

        top_filled_raster = raster_set[hu_raster_files[i]]['filled']
        bottom_filled_raster = raster_set[hu_raster_files[i + 1]]['filled']

        # Modify properties of layer using merge up logical array
        # loops over 'top' rasters for layers above current one.
        for j in xrange(0, i, 2):
            mod_j_idx = j / 2

            # Identify intersection of where this layer is inactive and where layer below is active
            # Modify top of mesh in this layer to the merge up locations,
            # but excluding those where the thickness mask was applied
            merge_up_locations = ~(cumulative_active[mod_j_idx] | inactive_cells[mod_i_idx])

            # Find where existing mesh is higher than candidate replacement from
            # next layer, so we know where to keep existing mesh
            new_elev = (top_filled_raster - mesh[mod_j_idx]) > 0.0

            # Find cells where they need updating and where cells aren't already
            # defined by zone, i.e. inactive cells
            modify_cells = new_elev & ~(cumulative_active[mod_j_idx] | ~inactive_cells[mod_i_idx])

            # Combine the two above logical arrays to apply in mesh modification
            combined = (merge_up_locations | modify_cells)
            if j > 0:
                # Check above layer which might have active cells for which we don't want
                # to modify the bottom of the layer
                combined = combined & ~cumulative_active[mod_j_idx - 1]
            # End if

            mesh[mod_j_idx][combined] = top_filled_raster[combined]

            zone_matrix[mod_j_idx][merge_up_locations] = mod_i_idx + 1  # Setting Zone ID to bottom raster id

            top_merge_locs = top_filled_raster[merge_up_locations]
            bot_merge_locs = bottom_filled_raster[merge_up_locations]
            merge_diff = top_merge_locs - bot_merge_locs
            next_j = mod_j_idx + 1
            proposed_heights = (top_merge_locs - (next_j * (merge_diff / next_j)))
            # proposed_heights = bot_merge_locs + ((top_merge_locs - bot_merge_locs) * 0.5)

            mesh[mod_j_idx + 1][merge_up_locations] = proposed_heights

            # After modifying the layer, update cumulative active for modified layer
            cumulative_active[mod_j_idx] = cumulative_active[mod_j_idx] | ~inactive_cells[mod_i_idx]
        # End for

        # Set properties for layer from i/2
        zone_matrix[mod_i_idx][inactive_cells[mod_i_idx]] = -1
        zone_matrix[mod_i_idx][~inactive_cells[mod_i_idx]] = mod_i_idx + 1
        # Set the base for layer i/2

        # Implement bottom of mesh for this layer BUT first check that the
        # inactive cell levels are not higher than the top of layer.
        # First identify where such cells might occur

        level_test = (bottom_filled_raster - mesh[mod_i_idx]) < 0

        # Get array of inactive cells in layer:
        # layer_inactive = ~cumulative_active[i/2]
        # Identify logical array of cells to modify:
        # modify = level_test & layer_inactive
        mesh[mod_i_idx + 1][level_test] = bottom_filled_raster[level_test]

    # End for

    # Use the fill function to clean up the inactive parts of the mesh to assist
    # with the rendering
    for idx, mod_idx, raster in top_rasters:
        invalid_vals = np.ma.masked_array(zone_matrix[mod_idx], zone_matrix[mod_idx] == -1).mask

        if idx == 0:
            mesh[idx] = fill(mesh[idx], invalid=np.ma.masked_array(zone_matrix[idx], zone_matrix[idx] == -1).mask)
        else:
            prev_invalid = np.ma.masked_array(zone_matrix[mod_idx - 1], zone_matrix[mod_idx - 1] == -1).mask
            mesh[mod_idx] = fill(mesh[mod_idx], invalid=invalid_vals | prev_invalid)
        # End if

        if mod_idx == (num_layers - 1):
            mesh[mod_idx + 1] = fill(mesh[mod_idx + 1], invalid=invalid_vals)
    # End for

    # Check that all layers are not intersecting
    for idx, mod_idx, raster in top_rasters:
        test = (mesh[mod_idx] - mesh[mod_idx + 1]) < min_height
        if np.any(test == True):
            print 'Error when processing layer', mod_idx

    # Calculate the thickness of each layer on the final mesh
    for idx, mod_idx, raster in top_rasters:
        thickness[mod_idx] = mesh[mod_idx] - mesh[mod_idx + 1]
        zone_matrix[mod_idx][thickness[mod_idx] <= 0] = -1
    # End for

    zone_matrix = reclassIsolatedCells(zone_matrix)

    tester2 = True
    if tester2:
        raster_thickness = {}

        for idx, mod_idx, raster in top_rasters:
            this_raster = raster_set[raster]['filled']
            next_raster = raster_set[hu_raster_files[idx + 1]]['filled']
            raster_thickness[mod_idx] = np.ma.masked_where((zone_matrix[mod_idx] != mod_idx + 1),
                                                           this_raster - next_raster)
            print 'Using', raster, ' and ', hu_raster_files[idx + 1]
        # End for

        mesh_zone_thickness = {}

        for i in xrange(num_hu_files / 2):
            top = np.zeros((len(raster_set[hu_raster_files[0]]['original']),
                            len(raster_set[hu_raster_files[0]]['original'][0])))
            bot = np.zeros((len(raster_set[hu_raster_files[0]]['original']),
                            len(raster_set[hu_raster_files[0]]['original'][0])))

            bot[zone_matrix[i] == i + 1] = mesh[i + 1][zone_matrix[i] == i + 1]

            for layer in reversed(range(zone_matrix.shape[0])):
                top[zone_matrix[layer] == i + 1] = mesh[layer][zone_matrix[layer] == i + 1]

            mesh_zone_thickness[i] = np.ma.masked_where((zone_matrix[i] != i + 1), top - bot)
        # End for

        for i in xrange(thickness.shape[0]):
            fig = plt.figure(figsize=(10, 8))
            fig.add_subplot(1, 3, 1, aspect='equal')
            plt.imshow(mesh_zone_thickness[i].data, interpolation='none')
            # print(i, mesh_zone_thickness[i].min(), mesh_zone_thickness[i].max())
            plt.title('Thickness in mesh: ' + hu_raster_files[i * 2])
            plt.colorbar()
            fig.add_subplot(1, 3, 2, aspect='equal')
            plt.imshow(raster_thickness[i], interpolation='none')
            # print(i, raster_thickness[i].min(), mesh_zone_thickness[i].max())
            plt.title('Thickness in raster')
            plt.colorbar()
            fig.add_subplot(1, 3, 3, aspect='equal')
            plt.imshow(mesh_zone_thickness[i] - raster_thickness[i], interpolation='none')
            plt.title('Difference in thickness')
            plt.colorbar()

            if np.any((mesh_zone_thickness[i] - raster_thickness[i]) < 0):
                print("Negative thickness detected in raster {}".format(i))
                # print(mesh_zone_thickness[i])
                # print(raster_thickness[i])
        # End for

    tester = True
    if tester:
        # print thickness
        thickness2 = {}
        thick_zero = {}

        thick_shape = thickness.shape[0]

        for i in xrange(thick_shape):
            thickness2[i] = np.ma.masked_where((zone_matrix[i] < 0), thickness[i])
            thick_zero[i] = np.ma.masked_where((zone_matrix[i] != 0), thickness[i])
        # End for

        for i in xrange(thick_shape):
            fig = plt.figure()
            fig.add_subplot(1, 3, 1, aspect='equal')
            plt.imshow(thickness2[i], interpolation='none')
            plt.title('Thickness in active cells: {}'.format(hu_raster_files[i * 2][0:5]))
            plt.colorbar()
            fig.add_subplot(1, 3, 2, aspect='equal')
            plt.imshow(thick_zero[i], interpolation='none')
            plt.title('Areas where thickness is 0')
            plt.colorbar()
            fig.add_subplot(1, 3, 3, aspect='equal')
            plt.imshow(zone_matrix[i], interpolation='none')
            plt.title('Zonal delineation')
            plt.colorbar()
        # End for

    if np.any(thickness < 0.):
        print 'issues'
    # End if

    grid_info = [ncol, nrow, delc, delr, x0, y0]
    array2Vtk.build_vtk_from_array(grid_info, np.fliplr(mesh), ["z_elev"], [np.fliplr(mesh)], [
                                   "zone", "thickness"],
                                   [np.fliplr(zone_matrix), np.fliplr(thickness)],
                                   out_path, vtk_out)

    return mesh, zone_matrix, grid_info
# End map_raster_array_to_mesh()


if __name__ == "__main__":
    bottom_ext = "_2b"
    top_ext = "_1t"

    # hu_raster_path = r"C:\Workspace\part0075\MDB modelling\VAF_v2.0_ESRI_GRID\ESRI_GRID\Preprocessed_data\Hydrogeological_Unit_Layers\\"

    hu_raster_path = r"D:\development\campaspe\groundwater\devbox\integrated\data\structured_model_grid_5000m"
    # hu_raster_files = ["qa_1t_bb", "qa_2b_bb", "utqa_1t_bb", "utqa_2b_bb",
    # "utaf_1t_bb", "utaf_2b_bb", "lta_1t_bb", "lta_2b_bb", "cps_1t_bb",
    # "cps_2b_bb"]

    hu_raster_files = ["qa_1t_model_grid.tif", "qa_2b_model_grid.tif",
                       "utb_1t_model_grid.tif", "utb_2b_model_grid.tif",
                       "utqa_1t_model_grid.tif", "utqa_2b_model_grid.tif",
                       "utam_1t_model_grid.tif", "utam_2b_model_grid.tif",
                       "utaf_1t_model_grid.tif", "utaf_2b_model_grid.tif",
                       "lta_1t_model_grid.tif", "lta_2b_model_grid.tif",
                       "bse_1t_model_grid.tif", "bse_2b.tif_model_grid.tif"
                       ]
    dx = "5000m"

    #hu_raster_path = hu_raster_path[:-2] + dx + hu_raster_path[-2:]
    # hu_raster_files = [x[:-4] + '_' + dx + x[-4:] for x in hu_raster_files]

    #      "utb_1t_model_grid.bil", "utb_2b_model_grid.bil",
    #      "qa_1t_model_grid.bil", "qa_2b_model_grid.bil",   , "bse_1t_model_grid.bil", "bse_2b.tif_model_grid.bil"

    # hu_layers = {"qa":["qa_1t_mb", "qa_2b_mb"], "utb":["utb_1t_mb", "utb_2b_mb"], "utqa":["utqa_1t_mb", "utqa_2b_mb"], "utaf":["utaf_1t_mb", "utaf_2b_mb"], "lta":["lta_1t_mb", "lta_2b_mb"], "bse":["bse_1t_mb", "bse_2b_mb.tif"]}
    #hu_layers_ordering = ["qa", "utb", "utqa", "utam", "utaf", "lta", "bse"]
    hu_layers_ordering = ["qa", "utb", "utqa", "utam", "utaf", "lta", "bse"]
    # hu_ext = "_mb"
    vtk_out = 'Campaspe_all_model_mesh'

    out_path = hu_raster_path

    min_height = 1.0  # Minimum thickness of cells
    max_height = 1000  # Maximum thickness of cells

    ret = map_raster_array_to_mesh(hu_raster_path, hu_raster_files,
                                   out_path, vtk_out, min_height, max_height)
