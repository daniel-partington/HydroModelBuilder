import os

import numpy as np


class MeshGenerator(object):
    """
    Mesh geneartion methods.
    """

    def __init__(self, Model_Int, GIS_Int):
        """
        Constructor for MeshGenerator.

        :param Model_Int: instance of ModelInterface, ModelInterface object
        :param GIS_Int: instance of GISInterface, GISInterface object
        """
        self.ModelInterface = Model_Int
        self.GISInterface = GIS_Int

        self.mesh_type = None
        self.centroid2mesh2Dindex = None
        self.centroid2mesh3Dindex = None

        self.model_boundary = None
        self.data_boundary = None
        self.boundary_poly_file = None
        self.boundary_data_file = None

        self.model_mesh = None
        self.model_mesh_centroids = None
        self.mesh2centroid2Dindex = None
        self.model_mesh3D = None
        self.model_mesh3D_centroids = None
        self.model_layers = None
        self.model_features = None
        self.polyline_mapped = {}
        self.polygons_mapped = {}
        self.points_mapped = {}
        self.pilot_points = {}
    # End __init__()

    def build_centroids_array3D(self, array_ordering):
        """
        1. Builds an array of cell centroids to be used in interpolating from
        other points arrays onto the cell centroids for the 3D mesh

        2. Creates a dictionary with centroids as key and lay row col as entries
        If key isn't found then it returns nearest centroid
        """
        (X, Y) = self.model_mesh_centroids
        (lays, rows, cols) = self.model_mesh3D[0].shape

        x = X[0]
        y = [y[0] for y in Y]

        X = np.asarray([X] * (lays - 1))
        Y = np.asarray([Y] * (lays - 1))

        self.model_mesh3D_centroids = (X, Y, self.model_mesh3D[0])
        self.centroid2mesh3Dindex = {}
        if self.mesh_type == 'structured' and array_ordering.array_order == 'UL_RowColumn':
            for lay, row, col in [(lay, row, col) for lay in xrange(lays - 1)
                                  for row in xrange(rows) for col in xrange(cols)]:
                self.centroid2mesh3Dindex[(
                    x[col],
                    y[row],
                    (self.model_mesh3D[0][lay][row][col] +
                     self.model_mesh3D[0][lay + 1][row][col]) / 2.
                )] = [lay, row, col]
            # End for
        # End if
        return self.model_mesh3D_centroids
    # End build_centroids_array3D()

    def build_centroids_array(self, gridHeight, array_ordering):
        """
        1. Builds an array of cell centroids to be used in interpolating from
        other points arrays onto the cell centroids

        2. Creates a dictionary with centroids as key and row col as entries
        If key isn't found then it returns nearest centroid

        :param gridHeight: float, height of grid cells
        """
        (xmin, xmax, ymin, ymax) = self.model_boundary[0:4]
        gridHeight = gridHeight

        rows = int((ymax - ymin) / gridHeight) + 1  # get rows
        cols = int((xmax - xmin) / gridHeight) + 1  # get columns

        # Nasty workaround until I find out why extent is wrong:
        (xmin, xmax, ymin, ymax) = self.model_mesh.GetLayer().GetExtent()
        (xmin2, xmax2, ymin2, ymax2) = (xmin, xmin + cols * gridHeight, ymax - rows * gridHeight, ymax)

        x = np.linspace(xmin2 + gridHeight / 2.0, xmax2 - gridHeight / 2.0, cols)
        # To put it in MODFLOW ordering ... need to invoke the array ordering here
        # to handle automatically
        y = np.linspace(ymax2 - gridHeight / 2.0, ymin2 + gridHeight / 2.0, rows)

        X, Y = np.meshgrid(x, y)
        self.model_mesh_centroids = (X, Y)
        self.centroid2mesh2Dindex = {}
        self.mesh2centroid2Dindex = {}
        if self.mesh_type == 'structured' and array_ordering.array_order == 'UL_RowColumn':
            for row, col in [(r, c) for r in xrange(rows) for c in xrange(cols)]:
                self.centroid2mesh2Dindex[(x[col], y[row])] = [row, col]
                self.mesh2centroid2Dindex[(row, col)] = [x[col], y[row]]
            # End for
        # End if

        return self.model_mesh_centroids
    # End build_centroids_array()

    def build_3D_mesh_from_rasters(self, out_data_folder_grid, array_ordering, raster_files, raster_path,
                                   minimum_thickness, maximum_thickness, force=False):
        """
        :param raster_files: list, raster files to build mesh from
        :param raster_path: str, path to raster files
        :param minimum_thickness: float, minimum thickness
        :param maximum_thickness: float, maximum thickness
        :param force: bool, if True, forces a rebuild of the mesh
        """
        p_j = os.path.join
        mesh_pth = p_j(out_data_folder_grid, 'model_mesh.npy')
        zone_pth = p_j(out_data_folder_grid, 'zone_matrix.npy')
        if os.path.isfile(mesh_pth) & os.path.isfile(zone_pth) & ~force:
            print 'Using previously generated mesh'
            self.model_mesh3D = self.ModelInterface.load_array(mesh_pth), self.ModelInterface.load_array(zone_pth)
        else:
            self.model_mesh3D = self.GISInterface.build_3D_mesh_from_rasters(
                raster_files, raster_path, minimum_thickness, maximum_thickness)
            self.ModelInterface.save_array(mesh_pth, self.model_mesh3D[0])
            self.ModelInterface.save_array(zone_pth, self.model_mesh3D[1])
        # End if

        self.build_centroids_array3D(array_ordering)  # Build 3D centroids array
    # End build_3D_mesh_from_rasters()

    def reclassIsolatedCells(self, passes=1, assimilate=False):
        """
        Function to remove cells that are surrounded by non-active cells in the horizontal plane

        e.g. if cells with positive integer is surrounded in above, below and to each side, then reassign to -1.

        :param passes: int, number of passes to do
        :param assimilate: bool, whether to assimilate or not if surrounded by four of the same.
        """
        mesh3D_1 = self.model_mesh3D[1]

        def most_common_oneliner(L):
            return max(g(sorted(L)), key=lambda(x, v): (len(list(v)), -L.index(x)))[0]
        # End most_common_oneliner()

        # Clean up idle cells:
        (lay, row, col) = mesh3D_1.shape
        for p, k, j, i in [(p_i, k_i, j_i, i_i) for p_i in xrange(passes)
                           for k_i in xrange(lay) for j_i in xrange(row) for i_i in xrange(col)]:
            cell_zone = mesh3D_1[k][j][i]
            if cell_zone == -1:
                continue
            # End if

            target_zone = mesh3D_1[k]
            # Assimilate cell if surrounded by four of the same
            if assimilate:
                north = ((j > 0) and (target_zone[j - 1][i] != cell_zone))
                south = ((j < row - 1) and (target_zone[j + 1][i] != cell_zone))
                east = ((i < col - 1) and (target_zone[j][i + 1] != cell_zone))
                west = ((i > 0) and (target_zone[j][i - 1] != cell_zone))
                north_east = ((j > 0) and (i < col - 1) and (target_zone[j - 1][i + 1] != cell_zone))
                south_east = ((j < row - 1) and (i < col - 1) and (target_zone[j + 1][i + 1] != cell_zone))
                north_west = ((j > 0) and (i > 0) and (target_zone[j - 1][i - 1] != cell_zone))
                south_west = ((j < row - 1) and (i > 0) and (target_zone[j + 1][i - 1] != cell_zone))
                if (north and south and
                        east and west and
                        north_east and south_east and
                        north_west and south_west):
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

                    most_common = most_common_oneliner(neighbours)
                    if most_common != -1:
                        target_zone[j][i] = most_common_oneliner(neighbours)
                    # End if
                # End if
            # End if

            # Check North, South, East, West zones
            # If any condition is true, then continue on
            # Otherwise, set the cell to -1
            if (((j > 0) and (target_zone[j - 1][i] != -1)) or       # North
                ((j < row - 1) and (target_zone[j + 1][i] != -1)) or  # South
                ((i < col - 1) and (target_zone[j][i + 1] != -1)) or  # East
                    ((i > 0) and (target_zone[j][i - 1] != -1))):         # West
                continue
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
        # End for
    # End reclassIsolatedCells()

    def _hacky_centroid_to_mesh(self, gridHeight, centroid, dist_min, points_dist, cache={}):
        """
        Nasty (and slow!) workaround due to minor mismatch in centroids from mesh and separate
        generation in this class. Perhaps better to actually define this array in fishnet when
        define_structured_mesh is called

        :param gridHeight: float,
        :param centroid:
        :param dist_min: float, minimum distance
        :param points_dist: method, function to use to calculate distance between points
        :param cache: dict, cache to store previously calculated distance
        """
        closest_key = cache.get((centroid, dist_min), None)
        if closest_key is None:
            half_grid_height = gridHeight / 2.0
            for key in self.centroid2mesh2Dindex:
                dist = points_dist(centroid, key)
                if dist < dist_min:
                    dist_min, closest_key = dist, key
                    if dist_min < half_grid_height:
                        break
                    # End if
                # End if
            # End for

            cache[(centroid, dist_min)] = self.centroid2mesh2Dindex[closest_key]
        # End if

        return self.centroid2mesh2Dindex[closest_key]
    # End _hacky_centroid_to_mesh()

    def map_points_to_3Dmesh(self, kdtree, points, identifier=None):
        '''
        Function to map points to the 3D mesh

        :params kdtree: method, function to apply
        :params points:
        :params identifier:

        Returns: A dict with keys of all points (or identifiers) with each
        corresponding entry containing the i, j and k reference of the nearest
        cell center to the given point
        '''
        centroid2mesh3Dindex = self.centroid2mesh3Dindex
        model_mesh_points = np.array(centroid2mesh3Dindex.keys())

        if type(points) == list:
            points = np.array(points)
        # End if

        closest = kdtree(model_mesh_points, points)

        point2mesh_map = {}
        for index, point in enumerate(points):
            if identifier[0]:
                point2mesh_map[identifier[index]] = centroid2mesh3Dindex[
                    tuple(model_mesh_points[closest[index]])]
            else:
                point2mesh_map[point] = centroid2mesh3Dindex[
                    tuple(model_mesh_points[closest[index]])]
            # End if
        # End for

        return point2mesh_map
    # End map_points_to_3Dmesh()

    def map_obs_loc2mesh3D(self, observations, kdtree, method='nearest', ignore=[-1]):
        """
        This is a function to map the obs locations to the nearest node in the
        mesh
        """
        if method == 'nearest':
            for key in observations.obs_group:
                if observations.obs_group[key]['domain'] == 'porous':
                    points = [list(x) for x in observations.obs_group[
                        key]['locations'].to_records(index=False)]
                    observations.obs_group[key]['mapped_observations'] = \
                        self.map_points_to_3Dmesh(kdtree,
                                                  points,
                                                  identifier=observations.obs_group[key]['locations'].index)

                    # Check that 'mapped_observations' are in active cells and if not then set
                    # the observation to inactive
                    for obs_loc in observations.obs_group[key]['mapped_observations'].keys():
                        [k, j, i] = observations.obs_group[key] \
                            ['mapped_observations'][obs_loc]
                        if self.model_mesh3D[1][k][j][i] in ignore:
                            observations.obs_group[key]['time_series'].loc[
                                observations.obs_group[key]['time_series']
                                ['name'] == obs_loc, 'active'] = False
                        else:
                            observations.obs_group[key]['time_series'].loc[
                                observations.obs_group[key]
                                ['time_series']['name'] == obs_loc, 'zone'] = \
                                "{}{}".format(key, int(self.model_mesh3D[1][k][j][i]))
                        # End if

                elif observations.obs_group[key]['domain'] == 'surface':
                    if observations.obs_group[key]['real']:
                        points = [list(x) for x in observations.obs_group[
                            key]['locations'].to_records(index=False)]
                        observations.obs_group[key]['mapped_observations'] = self.map_points_to_2Dmesh(
                            points,
                            identifier=observations.obs_group[key]['locations'].index)

                        # Check that 'mapped_observations' are in active cells and if not then set
                        # the observation to inactive
                        for obs_loc in observations.obs_group[key]['mapped_observations'].keys():
                            [j, i] = observations.obs_group[key]['mapped_observations'][obs_loc]
                            if self.model_mesh3D[1][0][j][i] in ignore:
                                observations.obs_group[key]['time_series'] \
                                    .loc[observations
                                         .obs_group[key]['time_series']['name'] ==
                                         obs_loc, 'active'] = False
                            else:
                                observations.obs_group[key]['time_series'] \
                                    .loc[observations
                                         .obs_group[key]['time_series']['name'] ==
                                         obs_loc, 'zone'] = \
                                    "{}{}".format(key, int(
                                        self.model_mesh3D[1][k][j][i]))
                            # End if
                        # End for
                    else:
                        observations.obs_group[key]['mapped_observations'] = {}
                        for index, loc in enumerate(observations.obs_group[key]['locations']):
                            observations.obs_group[key]['mapped_observations'][index] = loc
                        # End for
                    # End if
                elif observations.obs_group[key]['domain'] == 'stream':
                    if observations.obs_group[key]['real']:
                        observations.obs_group[key]['mapped_observations'] = {}
                        for index, loc in enumerate(observations.obs_group[key]['locations']):
                            observations.obs_group[key]['mapped_observations'][index] = loc
                        # End for
                    else:
                        observations.obs_group[key]['mapped_observations'] = {}
                        for index, loc in enumerate(observations.obs_group[key]['locations']):
                            observations.obs_group[key]['mapped_observations'][index] = loc
                        # End for
                    # End if
                # End if
                ts = observations.obs_group[key]['time_series']
                ts = ts[ts['active'] == True]
                observations.obs_group[key]['time_series'] = ts
            # End for
        # End if

        # return observations
    # End map_obs_loc2mesh3D()
