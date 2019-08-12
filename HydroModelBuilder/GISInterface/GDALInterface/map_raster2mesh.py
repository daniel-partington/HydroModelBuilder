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
from skimage import measure

from . import array2Vtk


# Step 1. Load raster layers top and bottom
def plot_map(array):
    plt.figure()
    plt.imshow(array, interpolation='none')

def remove_detached_clusters(zone_array, inactive=-1, ignore=[], explain=False, 
                             save_folder=''):
    '''
    Function to make use of image processing to detect detached clusters of HGU
    cells which should be discluded from the mesh.

    :param zone_array: 3D numpy int array of zone which each mesh cell belongs to
    :param ignore: list of zone values to ignore, default is empty list
    
    Steps in the process:
    1. create copy of zone array and modify all positive integer zone values to 
    the same value while allowing ignoring of some zones which can be set to -1
    
    2. use image processing to determine clusters by connectivity
    
    3. reclass all cells in secondary clusters as -1, i.e. inactive
    
    '''
    
    zone_array_mod = zone_array.copy()
    zone_array_new = zone_array.copy()
    for ig in ignore:
        zone_array_mod[zone_array_mod==ig] = inactive

    for num in range(1, int(np.max(zone_array)) + 1):
        zone_array_mod[zone_array_mod == num] = 1

    zone_array_labeled = measure.label(zone_array_mod, connectivity=1, background=inactive)
    
    for label in np.unique(zone_array_labeled):
        # For secondary clusters (label > 1) reassign cell zone to inactive
        if label > 1:
            zone_array_new[np.where(zone_array_labeled == label)] = inactive

    if explain:
        for layer in range(zone_array.shape[0]):                     
            fig = plt.figure(figsize=(20,8))
            #plt.title('Using image processing to detect anomalies\n in automatically generated mesh')
            ax1 = fig.add_subplot(1, 4, 1)
            ax1.imshow(zone_array[layer], interpolation='none')
            ax1.set_title('1. Generated mesh where colours represent\n hydrogeological units (HGU)\n blue and red are bad')
            ax2 = fig.add_subplot(1, 4, 2)
            ax2.imshow(zone_array_mod[layer], interpolation='none')
            ax2.set_title('2. Simplified mesh where we differentiate\n useful (red) and bad cells (blue)\n')
            ax3 = fig.add_subplot(1, 4, 3)
            ax3.imshow(zone_array_labeled[layer], interpolation='none')
            ax3.set_title('3. Image processing to identify connectivity\n and detached useful cells to reassign\n')
            ax3 = fig.add_subplot(1, 4, 4)
            ax3.imshow(zone_array_new[layer], interpolation='none')
            ax3.set_title('4. Improved mesh without detached HGUs\n\n')
            plt.subplots_adjust(top=0.87, bottom=0.1, left=0.04, right=0.98)
            plt.savefig(os.path.join(save_folder, "Image_processing_for_mesh_anomaly_detection{}.png".format(layer)), dpi=300)
                           
    return zone_array_new

def reclassIsolatedCells(mesh3D_1, passes=1, assimilate=False):
    """Function to remove cells that are surrounded by non-active cells in the horizontal plane
    e.g. if cells with positive integer is surrounded in above, below and to each side, then reassign to -1.
    :param mesh3D_1:
        :param passes:  (Default value = 1)
    :param assimilate: Default value = False)
    :param passes:  (Default value = 1)
    """

    # Clean up idle cells:
    (lay, row, col) = mesh3D_1.shape
    for p in range(passes):
        for k in range(lay):
            for j in range(row):
                for i in range(col):
                    cell_zone = mesh3D_1[k][j][i]
                    if cell_zone == -1:
                        continue
                    # End if

                    target_zone = mesh3D_1[k]
                    # Assimilate cell if surrounded by four of the same
                    if assimilate:
                        if (((j > 0) and (target_zone[j - 1][i] != cell_zone)) and       # North
                            # South
                            ((j < row - 1) and (target_zone[j + 1][i] != cell_zone)) and
                            ((i < col - 1) and (target_zone[j][i + 1] != cell_zone)) and  # East
                            ((i > 0) and (target_zone[j][i - 1] != cell_zone)) and       # West
                            # North-East
                            ((j > 0) and (i < col - 1) and (target_zone[j - 1][i + 1] != cell_zone)) and
                            # South-East
                            ((j < row - 1) and (i < col - 1) and (target_zone[j + 1][i + 1] != cell_zone)) and
                            # North-West
                            ((j > 0) and (i > 0) and (target_zone[j - 1][i - 1] != cell_zone)) and
                                ((j < row - 1) and (i > 0) and (target_zone[j + 1][i - 1] != cell_zone))):          # South-West

                            neighbours = []
                            if j > 0:
                                neighbours += [target_zone[j - 1][i]]
                            if j < row - 1:
                                neighbours += [target_zone[j + 1][i]]
                            if i < col - 1:
                                neighbours += [target_zone[j][i + 1]]
                            if i > 0:
                                neighbours += [target_zone[j][i - 1]]
                            # end if
                            from itertools import groupby as g

                            def most_common_oneliner(L):
                                """
                                :param L:
                                """

                                return max(g(sorted(L)), key=lambda x_v: (len(list(x_v[1])), -L.index(x_v[0])))[0]

                            most_common = most_common_oneliner(neighbours)
                            if most_common != -1:
                                target_zone[j][i] = most_common_oneliner(neighbours)

                    # Check North, South, East, West zones
                    # If any condition is true, then continue on
                    # Otherwise, set the cell to -1
                    if (((j > 0) and (target_zone[j - 1][i] != -1)) or       # North
                        ((j < row - 1) and (target_zone[j + 1][i] != -1)) or  # South
                        ((i < col - 1) and (target_zone[j][i + 1] != -1)) or  # East
                            ((i > 0) and (target_zone[j][i - 1] != -1))):         # West
                        continue
                    # End if

                    # End if

                    #  Check neighbours
                    # if k > 0:
                    #    cell_above_zone = MM.GW_build[name].model_mesh3D[1][k-1][j][i]
                    #    if cell_above_zone != -1: #cell_zone:
                    #        cell_active = True
                    # if k < lay - 1:
                    #    cell_below_zone = MM.GW_build[name].model_mesh3D[1][k+1][j][i]
                    #    if cell_below_zone != -1: # == cell_zone:
                    #        cell_active = True

                    # None of the above conditions were true
                    target_zone[j][i] = -1
                    print("reclassifying")

                # End for
            # End for
        # End for
    # End for
    return mesh3D_1
# End reclassIsolatedCells()

# This solution taken from:
# http://stackoverflow.com/questions/3662361/fill-in-missing-values-with-nearest-neighbour-in-python-numpy-masked-arrays


def fill(data, invalid=None):
    """Replace the value of invalid 'data' cells (indicated by 'invalid')
    by the value of the nearest valid data cell
    Input:
        data:    numpy array of any dimension
        invalid: a logical array of same shape as 'data'. True cells set
                 where data value should be replaced.
                 If None (default), use: invalid  = np.isnan(data)
    Output:
        Return a filled array.
    :param data:
    :param invalid:  (Default value = None)
    """

    if invalid is None:
        invalid = np.isnan(data)

    ind = nd.distance_transform_edt(invalid, return_distances=False, return_indices=True)
    return data[tuple(ind)]


def map_raster_array_to_mesh(hu_raster_path, hu_raster_files, out_path, vtk_out,
                             min_height, max_height, plot=False):
    """
    :param hu_raster_path:
    :param hu_raster_files:
    :param out_path:
    :param vtk_out:
    :param min_height:
    :param max_height:
    """
    raster_set = {}

    print(("HU Raster Path:", hu_raster_path))
    for raster in hu_raster_files:
        fname = os.path.join(hu_raster_path, raster)  # + hu_ext
        print('Processing: ', fname)
        ds = gdal.Open(fname, gdalconst.GA_ReadOnly)
        raster_set[raster] = [ds.GetRasterBand(1).ReadAsArray(),
                              ds.GetRasterBand(1).GetNoDataValue()]

        ds = gdal.Open(fname, gdalconst.GA_ReadOnly)
        raster_set[raster] = [ds.GetRasterBand(1).ReadAsArray(),
                              ds.GetRasterBand(1).GetNoDataValue()]

        if raster == hu_raster_files[0]:
            ncol = ds.RasterXSize
            nrow = ds.RasterYSize
            (x0, delr, z0, y0, z1, delc) = ds.GetGeoTransform()
            delc = -delc
        ds = None

    for raster in hu_raster_files:
        data = raster_set[raster][0]
        invalid = np.ma.masked_equal(data, np.min(data))
        raster_set[raster] += [invalid]
        a = fill(data, invalid=invalid.mask)
        raster_set[raster] += [a]

        # Clean_up
        data = None
        invalid = None
        a = None
    # End for

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

    for i in range(len(hu_raster_files)):
        if i % 2 == 1:
            continue
        cumulative_active[i / 2] = ~mask_active[i / 2]
    # end for

    for i in range(index + 2, len(hu_raster_files)):
        if i % 2 == 1:
            continue

        # Modify properties of layer using merge up logical array
        for j in range(0, i):
            if j % 2 == 1:
                continue

            # Identify intersection of where this layer is inactive and where
            # layer below is active
            merge_up_locations = ~(cumulative_active[j / 2] | mask_active[i / 2])
            # Modify top of mesh in this layer to the merge up locations,
            # but excluding those where the thickness mask was applied

            # Find where existing mesh is higher than candidate replacement from
            # next layer, so we know where to keep existing mesh
            new_elev = raster_set[hu_raster_files[i]][3] - mesh[j / 2] > 0.

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
                test = mesh[j / 2 - 1] - mesh[j / 2] < 0.
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

            test = mesh[j / 2] - mesh[j / 2 + 1] < 0.
            # if np.any(test == True):
            #    print 'Error when processing layer test3', i/2

            # After modifying the layer, update cumulative active for modified layer
            cumulative_active[j / 2] = cumulative_active[j / 2] | ~mask_active[i / 2]

        # Set properties for layer from i/2
        zone_matrix[i / 2][mask_active[i / 2]] = -1
        zone_matrix[i / 2][~mask_active[i / 2]] = i / 2 + 1
        # Set the base for layer i/2

        # Implement bottom of mesh for this layer BUT first check that the
        # inactive cell levels are not higher than the top of layer.
        # First identify where such cells might occur
        level_test = raster_set[hu_raster_files[i + 1]][3] - mesh[i / 2] < 0.
        # Get array of inactive cells in layer:
        # layer_inactive = ~cumulative_active[i/2]
        # Identify logical array of cells to modify:
        # modify = level_test & layer_inactive

        # mesh[i / 2 + 1][level_test] = raster_set[hu_raster_files[i + 1]][3][level_test]
        mesh[i / 2 + 1] = raster_set[hu_raster_files[i + 1]][3]

    # End for

    # Use the fill function to clean up the inactive parts of the mesh to assist
    # with the rendering
    cumulative_active_all = cumulative_active[0]
#    for lay in range(cumulative_active.shape[0]):
#        plot_map(cumulative_active[lay])

    for lay in range(mesh.shape[0]):
        if lay == 0:
            mesh[lay] = fill(mesh[lay], invalid=np.ma.masked_array(
                zone_matrix[lay], zone_matrix[lay] == -1).mask)
        elif lay < mesh.shape[0] - 1:
            mesh[lay] = fill(mesh[lay], invalid=np.ma.masked_array(
                zone_matrix[lay - 1], zone_matrix[lay - 1] == -1).mask)
        else:
            mesh[lay] = fill(mesh[lay], invalid=np.ma.masked_array(
                zone_matrix[lay - 1], zone_matrix[lay - 1] == -1).mask)

    # Check that all layers are not intersecting
    for index, raster in enumerate(hu_raster_files):
        if index % 2 == 1:
            continue
        test = mesh[index / 2] - mesh[index / 2 + 1] <= min_height
        if np.any(test == True):
            print('Minimum thickness not met when processing layer ', index / 2)

    # Calculate the thickness of each layer on the final mesh
    for index, raster in enumerate(hu_raster_files):
        # Only work on even index
        if index % 2 == 1:
            continue

        thickness[index / 2] = mesh[index / 2] - mesh[index / 2 + 1]
        if plot:
            plot_map(np.ma.masked_where(thickness[index / 2] <= 0, thickness[index / 2]))
        zone_matrix[index / 2][thickness[index / 2] <= 0] = -1
    # End for

    zone_matrix = reclassIsolatedCells(zone_matrix)

    zone_matrix = remove_detached_clusters(zone_matrix, ignore=[7])

    tester2 = True
    if tester2:
        raster_thickness = {}

        for index, raster in enumerate(hu_raster_files):
            if index % 2 == 1:
                continue
            raster_thickness[index / 2] = np.ma.masked_where(
                (zone_matrix[index / 2] != index / 2 + 1), raster_set[raster][3] -
                raster_set[hu_raster_files[index + 1]][3])
            print('Using', raster, ' and ', hu_raster_files[index + 1])

        mesh_zone_thickness = {}

        for i in range(len(hu_raster_files) / 2):
            #            top = np.zeros((len(raster_set[hu_raster_files[0]][0]),
            #                            len(raster_set[hu_raster_files[0]][0][0])))
            #            bot = np.zeros((len(raster_set[hu_raster_files[0]][0]),
            #                            len(raster_set[hu_raster_files[0]][0][0])))
            top = np.zeros_like(mesh[0])
            bot = np.zeros_like(mesh[0])
            bot[zone_matrix[i] == i + 1] = mesh[i + 1][zone_matrix[i] == i + 1]

            for layer in reversed(list(range(zone_matrix.shape[0]))):
                top[zone_matrix[layer] == i + 1] = mesh[layer][zone_matrix[layer] == i + 1]

            mesh_zone_thickness[i] = np.ma.masked_where((zone_matrix[i] != i + 1), top - bot)

        if plot:
            for i in range(thickness.shape[0]):
                fig = plt.figure()
                fig.add_subplot(1, 3, 1, aspect='equal')
                plt.imshow(mesh_zone_thickness[i], interpolation='none')
                print((i, mesh_zone_thickness[i].min(), mesh_zone_thickness[i].max(), mesh_zone_thickness[i].mean()))
                plt.title('Thickness in mesh: ' + hu_raster_files[i * 2])
                plt.colorbar()
                fig.add_subplot(1, 3, 2, aspect='equal')
                plt.imshow(raster_thickness[i], interpolation='none')
                print((i, raster_thickness[i].min(), raster_thickness[i].max(), raster_thickness[i].mean()))
                plt.title('Thickness in raster')
                plt.colorbar()
                fig.add_subplot(1, 3, 3, aspect='equal')
                plt.imshow(mesh_zone_thickness[i] - raster_thickness[i], interpolation='none')
                plt.title('Difference in thickness')
                plt.colorbar()
            # End for

    tester = True
    if tester:
        # print thickness
        thickness2 = {}
        thick_zero = {}

        thick_shape = thickness.shape[0]

        for i in range(thick_shape):
            thickness2[i] = np.ma.masked_where((zone_matrix[i] < 0), thickness[i])
            thick_zero[i] = np.ma.masked_where((thickness[i] != 0), thickness[i])
        # End for
        if plot:
            for i in range(thick_shape):
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
        print('issues')
        print(np.where(thickness < 0.))
    # End if
    if plot:
        for lay in range(thickness.shape[0]):
            plot_map(np.ma.masked_where(thickness[lay] > 0., thickness[lay]))

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

    hu_raster_path = r"C:\Workspace\part0075\MDB modelling\testbox\00_Campaspe_Cascade\01_steady_state\structured_model_grid_100m"
    hu_raster_files = ["qa_1t_model_grid.tif", "qa_2b_model_grid.tif",
                       "utb_1t_model_grid.tif", "utb_2b_model_grid.tif",
                       "utqa_1t_model_grid.tif", "utqa_2b_model_grid.tif",
                       "utam_1t_model_grid.tif", "utam_2b_model_grid.tif",
                       "utaf_1t_model_grid.tif", "utaf_2b_model_grid.tif",
                       "lta_1t_model_grid.tif", "lta_2b_model_grid.tif",
                       "bse_1t_model_grid.tif", "bse_2b.tif_model_grid.tif"
                       ]

#    hu_raster_path = r"C:\Workspace\part0075\MDB modelling\testbox\Campaspe_River_Small\steady_state\structured_model_grid_100m"
#    hu_raster_files = ["qa_1t_model_grid.tif", "qa_2b_model_grid.tif",
#                       "utb_1t_model_grid.tif", "utb_2b_model_grid.tif",
#                       "utqa_1t_model_grid.tif", "utqa_2b_model_grid.tif"]

#    dx = "100m"
    dx = "100m"

    vtk_out = 'Campaspe_all_model_mesh'

    out_path = hu_raster_path  # r"C:\Workspace\part0075\MDB modelling\integrated\Modules\Groundwater\model_files\testing\\"

    # Minimum thickness of cells
    min_height = 0.0001
    # Maximum thickness of cells
    max_height = 1000

    ret = map_raster_array_to_mesh(hu_raster_path, hu_raster_files,
                                   out_path, vtk_out, min_height, max_height)
