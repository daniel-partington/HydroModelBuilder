import cPickle as pickle
import os
import shutil
import sys

import numpy as np
import pandas as pd
from scipy import spatial

from Utilities import interpolation
from Utilities.PilotPoints import pilotpoints


class GWModelBuilder(object):
    """
    The ModelBuilder class contains a number of useful tools for building
    numerical groundwater models in packages such as MODFLOW, by dealing with
    spatial data using GIS type objects for reading and manipulating spatial
    data before mapping it to a model grid and converting it into an easily readable
    array structure to pass to the groundwater model.
    """

    def __init__(self, name=None, model_type=None, mesh_type=None,
                 units=None, data_folder=None, out_data_folder=None, model_data_folder=None,
                 GISInterface=None, data_format='binary', target_attr=None, **kwargs):
        """
        :param name: Model name.
        :param model_type: Type of model, e.g. MODFLOW (makes use of flopy), HGS (uses pyHGS ... which is not developed yet ;) ).
        :param mesh_type: type of mesh to be used in the model, i.e. structured or unstructured.
        :param data_folder: Path to get model data from.
        :param out_data_folder: Path to store processed model data.
        :param model_boundary: Boundary for the model domain.
        :param model_mesh: Mesh for the model.

        :param

        *** This list and every function description needs significant updating ...
        """
        # Define the constants for the model data types to use for checking input
        self.types = ModelBuilderType()

        try:
            # -- tests to alert user to incorrect inputs ...
            assert model_type in self.types.model_types, "Model types must be of type: {}".format(
                self.types.model_types)
            assert mesh_type in self.types.mesh_types, "'Mesh types must be of type: {}".format(
                self.types.mesh_types)
            assert data_format in self.types.data_formats, "Data format must be of type: {}".format(
                self.types.data_formats)

        except AssertionError as e:
            import traceback
            _, _, tb = sys.exc_info()
            # traceback.print_tb(tb) # Fixed format
            tb_info = traceback.extract_tb(tb)
            filename, line, func, text = tb_info[-1]
            sys.exit("An error occured in {} on line {} with the message '{}'".format(filename, line, e))
        # End try

        self.name = name
        self.model_type = model_type
        self.mesh_type = mesh_type

        self.units = {}
        if units != None:
            if type(units) == list:
                self.set_units(length=units[0], time=units[1], mass=units[2])
            elif type(units) == dict:
                self.set_units(length=units['length'], time=units['time'], mass=units['mass'])
            # end if
        else:
            self.set_units(length='m', time='d', mass='kg')
        # End if
        self.data_folder = data_folder
        self.out_data_folder = out_data_folder
        self.model_data_folder = model_data_folder
        self.GISInterface = GISInterface
        self.data_format = data_format

        self.array_ordering = ArrayOrdering()
        if self.model_type.lower() == 'modflow':
            self.array_ordering.SetModflowArrays()
        elif self.model_type.lower() == 'hgs':
            self.array_ordering.SetHGSArrays()
        else:
            pass
        # End if

        self.properties = ModelProperties()
        self.boundaries = ModelBoundaries()
        self.parameters = ModelParameters()
        self.observations = ModelObservations()
        self.initial_conditions = ModelInitialConditions()

        # OTHER variables
        self.out_data_folder_grid = None
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
        self.model_time = ModelTime()
        self.polyline_mapped = {}
        self.polygons_mapped = {}
        self.points_mapped = {}
        self.pilot_points = {}

        # Create registers to store details of processed and imported data for
        # quick future loading
        self.model_register = []
        self.base_data_register = []
        self.gridded_data_register = []

        # Set default target attributes
        if target_attr is None:
            self.target_attr = [
                'name',
                'model_type',
                'mesh_type',
                'units',
                'data_folder',
                'out_data_folder',
                'model_data_folder',
                'data_format',
                'array_ordering',
                'boundaries',
                'properties',
                'parameters',
                'observations',
                'initial_conditions',
                'out_data_folder_grid',
                'model_boundary',
                'boundary_poly_file',
                'boundary_data_file',
                'model_time',
                'model_mesh_centroids',
                'mesh2centroid2Dindex',
                'model_mesh3D',
                'model_mesh3D_centroids',
                'model_layers',
                'model_features',
                'polyline_mapped',
                'polygons_mapped',
                'points_mapped',
                'pilot_points',
                'model_register',
                'base_data_register',
                'gridded_data_register',
                # Some necessary parameters for now which should be replaced later
                'gridHeight',
                'gridWidth'
            ]
        else:
            self.target_attr = target_attr
        # End if

        # Set all other kwargs as class attributes
        for key, value in kwargs.items():
            setattr(self, key, value)
        # End For

        # Set all key word arguments as attributes
        if self.GISInterface is not None:
            for key, value in self.__dict__.items():
                if type(value) is not object:
                    setattr(self.GISInterface, key, value)
                # End if
            # End for
        # end if

    def updateGISinterface(self):
        for key, value in self.__dict__.items():
            if type(value) is not object:
                setattr(self.GISInterface, key, value)

    def set_units(self, length='m', time='s', mass='kg'):

        if length in self.types.length:
            self.units['length'] = length
        else:
            print 'Length unit not recognised, please use one of: ', self.types.length
            print 'Default unit "m" is set'
        # end if
        if time in self.types.time:
            self.units['time'] = time
        else:
            print 'Time unit not recognised, please use one of: ', self.types.time
            print 'Default unit "s" is set'
        # end if
        if mass in self.types.mass:
            self.units['mass'] = mass
        else:
            print 'Mass unit not recognised, please use one of: ', self.types.mass
            print 'Default unit "kg is set'
        # end if

    def check_for_existing(self, f):
        """
        Function to determine if input files have previously been processed
        and if so to do nothing unless flagged otherwise. This is done by
        checking the output data path.
        """
        filename_suffixes = ['_model', '_grid']

        for suffix in filename_suffixes:
            if os.path.isfile(os.path.join(self.out_data_folder, f[:-4] +
                                           suffix + f[-4:])):
                print 'found processed file'
        # end for
        # if any(os.path.isfile(x) in filename for x in filename_suffixes):
        #    pass

    def save_array(self, filename, array):
        if self.data_format == self.types.data_formats[0]:  # if it is 'ascii'
            np.savetxt(filename, array)
        elif self.data_format == self.types.data_formats[1]:  # if it is 'binary'
            np.save(filename, array)
        else:
            print 'Data format not recognised, use "binary" or "ascii"'
        # end if

    def load_array(self, array_file):
        if array_file.endswith('txt'):
            return np.loadtxt(array_file)
        elif array_file.endswith('npy') or array_file.endswith('npz'):
            return np.load(array_file)
        else:
            print 'File type not recognised as "txt", "npy" or "npz" \n'
            sys.exit(1)
        # end if

    def save_dataframe(self, filename, df):
        df.to_hdf(filename + '.h5', 'table')

    def load_dataframe(self, filename):
        if filename.endswith('.h5'):
            return pd.read_hdf(filename, 'table')
        else:
            print 'File type not recognised as "h5"'
            sys.exit(1)
        # end if

    def save_obj(self, obj, filename):
        with open(filename + '.pkl', 'wb') as f:
            pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

    def load_obj(self, filename):
        if filename.endswith('.pkl'):
            with open(filename, 'rb') as f:
                print f
                print filename
                p = pickle.load(f)
                return p

        else:
            print 'File type not recognised as "pkl"'
            sys.exit(1)
        # end if

    def flush(self, mode=None):
        if mode == 'data':
            folder = self.out_data_folder
        elif mode == 'model':
            folder = self.model_data_folder
        else:
            print 'Expected mode to be either "data" or "model" but got: ', mode
            sys.exit(1)

        if folder == None:
            print 'No folder set, so no flushing'
            sys.exit(1)
        # end if
        for the_file in os.listdir(folder):
            file_path = os.path.join(folder, the_file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
                # end if
            except Exception as e:
                print(e)
        # End for

    def set_model_boundary_from_corners(self, xmin, xmax, ymin, ymax):
        """
        Function to set model boundary based on x and y bounds: xmin, xmax, ymin, ymax

                    ._________.(xmax, ymax)
                    |         |
                    |         |
                    |         |
                    ._________.
         (xmin,ymin)
        """

        # self.model_boundary = [xmin, xmax, ymin, ymax]
        self.model_boundary[0] = xmin
        self.model_boundary[1] = xmax
        self.model_boundary[2] = ymin
        self.model_boundary[3] = ymax

    def set_model_boundary_from_polygon_shapefile(self, shapefile_name, shapefile_path=None):
        if shapefile_path is None:
            shapefile_path = self.data_folder
        # end if

        self.model_boundary, self.boundary_poly_file = self.GISInterface.set_model_boundary_from_polygon_shapefile(
            shapefile_name, shapefile_path)
        return self.model_boundary, self.boundary_poly_file

    def set_data_boundary_from_polygon_shapefile(self, shapefile_name=None, shapefile_path=None, buffer_dist=None):
        if shapefile_name is None:
            shapefile_name = self.boundary_poly_file
        if shapefile_path is None:
            shapefile_path = self.out_data_folder
        # end if
        self.data_boundary, self.boundary_data_file = self.GISInterface.set_data_boundary_from_polygon_shapefile(
            shapefile_name, shapefile_path=shapefile_path, buffer_dist=buffer_dist)
        self.updateGISinterface()
        return self.data_boundary, self.boundary_data_file

    def build_centroids_array(self, gridHeight):
        """
        1. Builds an array of cell centroids to be used in interpolating from
        other points arrays onto the cell centroids

        2. Creates a dictionary with centroids as key and row col as entries
        If key isn't found then it returns nearest centroid
        """

        (xmin, xmax, ymin, ymax) = self.model_boundary[0:4]
        gridHeight = gridHeight

        # get rows
        rows = int((ymax - ymin) / gridHeight) + 1
        # get columns
        cols = int((xmax - xmin) / gridHeight) + 1

        # Nasty workaround until I find out why extent is wrong:
        (xmin, xmax, ymin, ymax) = self.model_mesh.GetLayer().GetExtent()
        (xmin2, xmax2, ymin2, ymax2) = (xmin, xmin + cols * gridHeight, ymax - rows * gridHeight, ymax)

        x = np.linspace(xmin2 + gridHeight / 2.0, xmax2 - gridHeight / 2.0, cols)
        # y = np.linspace(ymin+gridHeight/2.0, ymax-gridHeight/2.0, rows)
        # To put it in MODFLOW ordering ... need to invoke the array ordering here
        # to handle automatically
        y = np.linspace(ymax2 - gridHeight / 2.0, ymin2 + gridHeight / 2.0, rows)

        X, Y = np.meshgrid(x, y)
        self.model_mesh_centroids = (X, Y)

        self.centroid2mesh2Dindex = {}
        self.mesh2centroid2Dindex = {}

        if self.mesh_type == 'structured':
            if self.array_ordering.array_order == 'UL_RowColumn':
                for row in range(rows):
                    for col in range(cols):
                        self.centroid2mesh2Dindex[(x[col], y[row])] = [row, col]
                        self.mesh2centroid2Dindex[(row, col)] = [x[col], y[row]]
                    # end for
                # end for
            # end if
        # end if

        return self.model_mesh_centroids

    def build_centroids_array3D(self):
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

        X = np.asarray([X] * (lays - 1))  # np.repeat(X[np.newaxis, :, :], (lays-1), axis=0)
        Y = np.asarray([Y] * (lays - 1))  # np.repeat(Y[np.newaxis, :, :], (lays-1), axis=0)

        self.model_mesh3D_centroids = (X, Y, self.model_mesh3D[0])

        self.centroid2mesh3Dindex = {}

        if self.mesh_type == 'structured':
            if self.array_ordering.array_order == 'UL_RowColumn':
                for lay in xrange(lays - 1):
                    for row in xrange(rows):
                        for col in xrange(cols):
                            self.centroid2mesh3Dindex[(
                                x[col],
                                y[row],
                                (self.model_mesh3D[0][lay][row][col] +
                                 self.model_mesh3D[0][lay + 1][row][col]) / 2.
                            )] = [lay, row, col]
                        # end for
                    # end for
                # end for
            # end if
        # end if
        return self.model_mesh3D_centroids

    def define_structured_mesh(self, gridHeight, gridWidth):

        self.gridHeight = gridHeight
        self.gridWidth = gridWidth
        self.out_data_folder_grid = os.path.join(self.model_data_folder,
                                                 'structured_model_grid_{}m'.format(gridHeight))
        self.updateGISinterface()
        self.model_mesh = self.GISInterface.define_structured_mesh(gridHeight, gridWidth)
        self.model_mesh_centroids = self.build_centroids_array(self.gridHeight)
        return (self.out_data_folder_grid, self.gridHeight, self.gridWidth,
                self.model_mesh_centroids, self.model_mesh)

    def read_rasters(self, files, path=None):

        return self.GISInterface.read_rasters(files, path=path)

    def map_rasters_to_grid(self, raster_files, raster_path):

        return self.GISInterface.map_rasters_to_grid(raster_files, raster_path)

    def map_heads_to_mesh_by_layer(self):

        return self.GISInterface.map_heads_to_mesh_by_layer()

    def create_basement_bottom(self, hu_raster_path, surface_raster_file, basement_top_raster_file, basement_bot_raster_file, output_path, raster_driver='GTiff'):
        """
        Utility to build a bottom basement array where it doesn't exist based on top of bedrock, surface elevation and a thickness function

        writes a raster file for the bottom basement array

        ** This perhaps would better sit in a separate utilities folder ...
        """
        return self.GISInterface.create_basement_bottom(hu_raster_path, surface_raster_file,
                                                        basement_top_raster_file,
                                                        basement_bot_raster_file,
                                                        output_path,
                                                        raster_driver=raster_driver)

    def build_3D_mesh_from_rasters(self, raster_files, raster_path, minimum_thickness, maximum_thickness, force=False):
        p_j = os.path.join
        mesh_pth = p_j(self.out_data_folder_grid, 'model_mesh.npy')
        zone_pth = p_j(self.out_data_folder_grid, 'zone_matrix.npy')
        if os.path.isfile(mesh_pth) & os.path.isfile(zone_pth) & ~force:
            print 'Using previously generated mesh'
            self.model_mesh3D = self.load_array(
                mesh_pth), self.load_array(zone_pth)
        else:
            self.model_mesh3D = self.GISInterface.build_3D_mesh_from_rasters(
                raster_files, raster_path, minimum_thickness, maximum_thickness)
            self.save_array(mesh_pth, self.model_mesh3D[0])
            self.save_array(zone_pth, self.model_mesh3D[1])
        # end if

        # Build 3D centroids array:
        self.build_centroids_array3D()

    def reclassIsolatedCells(self, passes=1, assimilate=False):
        """
        Function to remove cells that are surrounded by non-active cells in the horizontal plane

        e.g. if cells with positive integer is surrounded in above, below and to each side, then reassign to -1.

        """

        # Clean up idle cells:
        (lay, row, col) = self.model_mesh3D[1].shape
        for p in xrange(passes):
            for k in xrange(lay):
                for j in xrange(row):
                    for i in xrange(col):
                        cell_zone = self.model_mesh3D[1][k][j][i]
                        if cell_zone == -1:
                            continue
                        # End if

                        target_zone = self.model_mesh3D[1][k]
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
                                    return max(g(sorted(L)), key=lambda(x, v): (len(list(v)), -L.index(x)))[0]

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

                    # End for
                # End for
            # End for
        # End for

    def read_poly(self, filename, path=None, poly_type='polyline'):

        return self.GISInterface.read_poly(filename, path, poly_type=poly_type)

    def map_polyline_to_grid(self, polyline_obj):

        if type(polyline_obj) is str:
            polyline_obj = self.read_poly(polyline_obj)
        # end if
        if os.path.sep in polyline_obj.GetDescription():
            poly_name = polyline_obj.GetDescription().split(os.path.sep)[-1]
        else:
            poly_name = polyline_obj.GetDescription()
        # end if

        if os.path.exists(os.path.join(self.out_data_folder_grid, poly_name + '_mapped.pkl')):
            print "Using previously mapped polyline to grid object"
            self.polyline_mapped[poly_name] = self.load_obj(os.path.join(self.out_data_folder_grid,
                                                                         poly_name + '_mapped.pkl'))
        else:
            temp = []
            self.polyline_mapped[poly_name] = self.GISInterface.map_polyline_to_grid(polyline_obj)
            for item in self.polyline_mapped[poly_name]:
                centroid = (float(item[1][0]), float(item[1][1]))
                # replace centroid with row col
                try:
                    grid_loc = self.centroid2mesh2Dindex[centroid]
                except:
                    # This is a very nasty workaround AND slow due to mismatch in
                    # centroids from mesh and separate generation in this class
                    # Perhaps better to actually define this array in fishnet
                    # when define structured mesh is called
                    def pointsdist(p1, p2):
                        return ((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)**0.5

                    dist_min = 1E6
                    for key in self.centroid2mesh2Dindex:
                        dist = pointsdist(centroid, key)
                        if dist < dist_min:
                            dist_min = dist
                            closest_key = key
                            if dist_min < self.gridHeight / 2.0:
                                break

                    grid_loc = self.centroid2mesh2Dindex[closest_key]

                temp += [[grid_loc, item[0]]]
            # end for
            self.polyline_mapped[poly_name] = temp
            self.save_obj(temp, os.path.join(self.out_data_folder_grid, poly_name + '_mapped'))
        # end if
        self.gridded_data_register += [poly_name]

    def map_polygon_to_grid(self, polygon_obj, feature_name=None):

        if type(polygon_obj) is str:
            polygon_obj = self.read_poly(polygon_obj)
        # end if
        if os.path.sep in polygon_obj.GetDescription():
            poly_name = polygon_obj.GetDescription().split(os.path.sep)[-1]
        else:
            poly_name = polygon_obj.GetDescription()
        # end if

        if os.path.exists(os.path.join(self.out_data_folder_grid, poly_name + '_mapped.pkl')):
            print "Using previously mapped polygons to grid object"
            self.polygons_mapped[poly_name] = self.load_obj(os.path.join(self.out_data_folder_grid,
                                                                         poly_name + '_mapped.pkl'))
        else: 
            self.polygons_mapped[poly_name] = \
                self.GISInterface.map_polygon_to_grid(polygon_obj, 
                                                      out_fname=os.path.join( \
                                                          self.out_data_folder_grid, poly_name + '_mapped'), 
                                                      pixel_size=self.model_mesh_centroids[0][0][1] - self.model_mesh_centroids[0][0][0], 
                                                      bounds=self.model_mesh.GetLayer().GetExtent(), 
                                                      feature_name=feature_name)
            
            self.save_obj(self.polygons_mapped[poly_name], 
                          os.path.join(self.out_data_folder_grid, 
                                       poly_name + '_mapped'))
        # end if
        
        self.gridded_data_register += [poly_name]
       

    def read_points_data(self, filename, path=None):

        return self.GISInterface.read_points_data(filename, path)

    def map_points_to_grid(self, points_obj, feature_id=None):

        if type(points_obj) is str:
            points_obj = self.read_points_data(points_obj)
        # end if
        if os.path.sep in points_obj.GetDescription():
            point_name = points_obj.GetDescription().split(os.path.sep)[-1]
        else:
            point_name = points_obj.GetDescription()
        # end if

        if os.path.exists(os.path.join(self.out_data_folder_grid, point_name + '_mapped.pkl')):
            print "Using previously mapped points to grid object"
            self.points_mapped[point_name] = self.load_obj(os.path.join(self.out_data_folder_grid,
                                                                        point_name + '_mapped.pkl'))
        else:
            temp = []
            self.points_mapped[point_name] = self.GISInterface.map_points_to_grid(
                points_obj, feature_id=feature_id)
            # Go from centroids to ij indexing
            for item in self.points_mapped[point_name]:
                centroid = (float(item[1][0]), float(item[1][1]))
                # replace centroid with row col
                try:
                    grid_loc = self.centroid2mesh2Dindex[centroid]
                except:
                    # This is a very nasty workaround AND slow due to minor mismatch in
                    # centroids from mesh and separate generation in this class
                    # Perhaps better to actually define this array in fishnet
                    # when define_structured_mesh is called
                    def pointsdist(p1, p2):
                        return ((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)**0.5

                    dist_min = 1E6
                    for key in self.centroid2mesh2Dindex:
                        dist = pointsdist(centroid, key)
                        if dist < dist_min:
                            dist_min = dist
                            closest_key = key
                            if dist_min < self.gridHeight / 2.0:
                                break
                    # print 'Closest key is: ', closest_key
                    grid_loc = self.centroid2mesh2Dindex[closest_key]

                temp += [[grid_loc, item[0]]]

            self.points_mapped[point_name] = temp
            self.save_obj(temp, os.path.join(self.out_data_folder_grid, point_name + '_mapped'))

        self.gridded_data_register += [point_name]

    def do_kdtree(self, model_mesh_points, points):
        mytree = spatial.cKDTree(model_mesh_points)
        dist, indexes = mytree.query(points)
        return indexes


    def map_points_to_2Dmesh(self, points, identifier=None):
        '''
        Function to map points to the 2D horizontal mesh (i.e. on the xy plane)

        Returns: A dict with keys of all points (or identifiers) with each
        corresponding entry containing the i and j reference of the nearest
        cell center to the given point
        '''
        model_mesh_points = np.array(self.centroid2mesh2Dindex.keys())

        if type(points) == list:
            points = np.array(points)
        # end if


        closest = self.do_kdtree(model_mesh_points, points)
        point2mesh_map = {}

        for index, point in enumerate(points):
            if identifier[0]:
                point2mesh_map[identifier[index]] = self.centroid2mesh2Dindex[
                    tuple(model_mesh_points[closest[index]])]
            else:
                point2mesh_map[point] = self.centroid2mesh2Dindex[
                    tuple(model_mesh_points[closest[index]])]
            # end if
        # end for

        return point2mesh_map

    def find_closest_points_between_two_lists(self, list1, list2):
        np_arr1 = np.array(list1)
        np_arr2 = np.array(list2)
        closest = self.do_kdtree(np_arr1, np_arr2)
        return closest

    def map_points_to_3Dmesh(self, points, identifier=None):
        '''
        Function to map points to the 3D mesh

        Returns: A dict with keys of all points (or identifiers) with each
        corresponding entry containing the i, j and k reference of the nearest
        cell center to the given point
        '''

        model_mesh_points = np.array(self.centroid2mesh3Dindex.keys())

        if type(points) == list:
            points = np.array(points)
        # end if

        def do_kdtree(model_mesh_points, points):
            mytree = spatial.cKDTree(model_mesh_points)
            dist, indexes = mytree.query(points)
            return indexes

        closest = do_kdtree(model_mesh_points, points)
        point2mesh_map = {}

        for index, point in enumerate(points):
            if identifier[0]:
                point2mesh_map[identifier[index]] = self.centroid2mesh3Dindex[
                    tuple(model_mesh_points[closest[index]])]
            else:
                point2mesh_map[point] = self.centroid2mesh3Dindex[
                    tuple(model_mesh_points[closest[index]])]
            # end if
        # end for

        return point2mesh_map

    def map_obs_loc2mesh3D(self, method='nearest', ignore=[-1]):
        """
        This is a function to map the obs locations to the nearest node in the
        mesh

        """

        if method == 'nearest':
            for key in self.observations.obs_group:
                if self.observations.obs_group[key]['domain'] == 'porous':
                    points = [list(x) for x in self.observations.obs_group[
                        key]['locations'].to_records(index=False)]
                    self.observations.obs_group[key]['mapped_observations'] = self.map_points_to_3Dmesh(
                        points, identifier=self.observations.obs_group[key]['locations'].index)

                    # Check that 'mapped_observations' are in active cells and if not then set
                    # the observation to inactive
                    for obs_loc in self.observations.obs_group[key]['mapped_observations'].keys():
                        [k, j, i] = self.observations.obs_group[key]['mapped_observations'][obs_loc]
                        if self.model_mesh3D[1][k][j][i] in ignore:
                            self.observations.obs_group[key]['time_series'].loc[
                                self.observations.obs_group[key]['time_series']['name'] == obs_loc, 'active'] = False
                        else:
                            self.observations.obs_group[key]['time_series'].loc[
                                self.observations.obs_group[key]['time_series']['name'] == obs_loc, 'zone'] = "{}{}".format(key, int(self.model_mesh3D[1][k][j][i]))

                elif self.observations.obs_group[key]['domain'] == 'surface':
                    if self.observations.obs_group[key]['real']:
                        points = [list(x) for x in self.observations.obs_group[
                            key]['locations'].to_records(index=False)]
                        self.observations.obs_group[key]['mapped_observations'] = self.map_points_to_2Dmesh(
                            points, identifier=self.observations.obs_group[key]['locations'].index)

                        # Check that 'mapped_observations' are in active cells and if not then set
                        # the observation to inactive
                        for obs_loc in self.observations.obs_group[key]['mapped_observations'].keys():
                            [j, i] = self.observations.obs_group[
                                key]['mapped_observations'][obs_loc]
                            if self.model_mesh3D[1][0][j][i] in ignore:
                                self.observations.obs_group[key]['time_series'].loc[
                                    self.observations.obs_group[key]['time_series']['name'] == obs_loc, 'active'] = False
                            else:
                                self.observations.obs_group[key]['time_series'].loc[
                                    self.observations.obs_group[key]['time_series']['name'] == obs_loc, 'zone'] = "{}{}".format(key, int(self.model_mesh3D[1][k][j][i]))
                            # end if
                        # end for
                    else:
                        self.observations.obs_group[key]['mapped_observations'] = {}
                        for index, loc in enumerate(self.observations.obs_group[key]['locations']):
                            self.observations.obs_group[key]['mapped_observations'][index] = loc
                        # end for
                    # end if
                # end if
                ts = self.observations.obs_group[key]['time_series']
                ts = ts[ts['active'] == True]
                self.observations.obs_group[key]['time_series'] = ts
            # end for
        # end if

    def getXYpairs(self, points_obj, feature_id=None):

        return self.GISInterface.getXYpairs(points_obj, feature_id=feature_id)

    def points_shapefile_obj2dataframe(self, points_obj, feature_id=None):
        """
        Converts GIS object to dataframe, which is designed for point style objects only
        """
        points_dict = self.GISInterface.getXYpairs(points_obj, feature_id=feature_id)
        df_points = pd.DataFrame.from_dict(points_dict, orient='index')
        df_points.columns = ['Easting', 'Northing']
        df_points[feature_id] = df_points.index
        return df_points

    def point_values_from_raster(self, points, raster_obj):
        """
        Get values from raster at given points defined by list or by points object
        points must be a list made up of point lists or points shapefile object
        """

        return self.GISInterface.point_values_from_raster(points, raster_obj)

    def map_points_to_raster_layers(self, points, depths, rasters):

        # create boolean array
        len_list = len(points)
        layers = len(rasters) / 2  # Intentional integer division
        points_layer = np.full((layers, len_list), False, dtype=bool)

        return self.GISInterface.map_points_to_raster_layers(points, depths, rasters, points_layer)

    def interpolate_points2mesh(self, points_obj, values_dataframe, feature_id=None,
                                method='nearest', use='griddata',
                                function='multiquadric', epsilon=2):

        if isinstance(points_obj, list) or isinstance(points_obj, np.ndarray):
            points = points_obj
            values = values_dataframe
        else:
            points_dict = self.getXYpairs(points_obj, feature_id=feature_id)
            points = []
            values = []
            for key in points_dict:
                points += [points_dict[key]]
                values += [float(values_dataframe[key])]

        return interpolation.Interpolator(self.mesh_type, np.array(points), np.array(values),
                                          self.model_mesh_centroids, method=method, use=use,
                                          function=function, epsilon=epsilon)

    def interpolate_points2meshByLayer():
        pass

    @staticmethod
    def _findInterval(row, times):
        key_time = row['datetime']
        lower_time = times[0]
        for period, time in enumerate(times):
            if period > 0:
                if lower_time <= key_time < time:
                    return period - 1
            lower_time = time
        return np.nan

    def map_obs2model_times(self):
        """
        This is a function to map the obs at different times within the bounding interval in the
        model times intervals

        """

        if self.model_time.t['steady_state']:
            for key in self.observations.obs_group.keys():
                self.observations.obs_group[key]['time_series']['interval'] = 0
        else:
            for key in self.observations.obs_group.keys():
                self.observations.obs_group[key]['time_series']['interval'] = self.observations.obs_group[key][
                    'time_series'].apply(lambda row: self._findInterval(row, self.model_time.t['dateindex']), axis=1)
                # remove np.nan values from the obs as they are not relevant
                if self.observations.obs_group[key]['real']:
                    self.observations.obs_group[key]['time_series'] = self.observations.obs_group[key][
                        'time_series'][pd.notnull(self.observations.obs_group[key]['time_series']['interval'])]
                    self.observations.obs_group[key]['time_series'] = self.observations.obs_group[key][
                        'time_series'][self.observations.obs_group[key]['time_series']['value'] != 0.]

    def updateModelParameters(self, fname, verbose=True):
        with open(fname, 'r') as f:
            text = f.readlines()
            # Remove header
            text = text[1:]
            # Read in parameters and replace values in parameters class for param
            updated = {}
            for key in self.parameters.param.keys():
                updated[key] = False

            for line in text:
                param_name, value = line.strip('\n').split('\t')
                value = value.lstrip()
                param_name = param_name.strip()
                if param_name in self.parameters.param.keys():
                    self.parameters.param[param_name]['PARVAL1'] = float(value)
                    updated[param_name] = True
                else:
                    if verbose:
                        print 'Parameter not defined in model: ', param_name

        if verbose:
            were_updated = [key for key in updated.keys() if updated[key] == True]
            if len(were_updated) > 0:
                print 'Parameters updated for : ', were_updated
            not_updated = [key for key in updated.keys() if updated[key] == False]
            if len(not_updated) > 0:
                print 'Parameters unchanged for : ', not_updated

    def create_pilot_points(self, name):
        self.pilot_points[name] = pilotpoints.PilotPoints(
            output_directory=self.out_data_folder_grid)

    def save_pilot_points(self):
        self.save_obj(self.pilot_points, os.path.join(self.out_data_folder_grid, 'pilot_points'))

    def load_pilot_points(self, fname):
        pp = self.load_obj(fname)
        for key in pp.keys():
            self.pilot_points[key] = pp[key]

    def add2register(self, addition):

        self.model_register += addition

    def writeRegister2file(self):

        with open(os.path.join(self.out_data_folder, 'model_register.dat')) as f:
            for item in self.model_register:
                f.write(item)

    def points2shapefile(self, points_array, shapefile_name):
        """ Needs writing ... """
        self.GISInterface.points2shapefile(points_array, shapefile_name)

    def mesh3DToVtk(self, val_array, val_name, out_path, vtk_out):
        '''
        Function to write the mesh array
        '''
        nrow, ncol = self.model_mesh3D[1][0].shape
        delc, delr = self.gridWidth, self.gridHeight
        x0, y0 = self.model_boundary[0], self.model_boundary[3]
        grid_info = [ncol, nrow, delc, delr, x0, y0]
        mesh, zone_matrix = self.mesh_array[0], self.mesh_array[1]
        from HydroModelBuilder.GISInterface.GDALInterface import array2Vtk

        array2Vtk.build_vtk_from_array(grid_info, np.fliplr(mesh), ["z_elev"],
                                       [np.fliplr(mesh)], ["zone", val_name],
                                       [np.fliplr(zone_matrix), np.fliplr(val_array)],
                                       out_path, vtk_out)

    def get_uppermost_active(self):
        # Not required at the moment
        pass

    def package_data(self):
        """
        Option to save all important attributes of GWModelBuilder class to
        allow quick loading of data that may have required transforms and
        processing in its orignial state.

        This will include GIS type objects
        """
        pass

    def package_model(self):
        """
        Option to save all important attributes of model class to allow quick
        loading and manipulation of model without separate data and grid generation
        building scripts

        This will not include any GIS type objects
        """

        target_attr = self.target_attr
        packaged_model = {k: self.__dict__[k] for k in self.__dict__ if k in target_attr}

        # Hack to fix up model boundary which contains a gdal object as well:
        packaged_model['model_boundary'] = packaged_model['model_boundary'][0:4]

        self.save_obj(packaged_model, os.path.join(self.out_data_folder_grid,
                                                   self.name + '_packaged'))

    ###########################################################################

    def __exit__(self):
        self.writeRegister2file()
        self.save_mapped_dictionaries()


class GeneralModelDataStructure(object):
    """
    This class describes a general model structure to be used for
    setting up the groundwater models. The class is made up of numpy
    arrays that describe the properties and forcing to be applied to the model

    """

    def __init__(self, node_type=None, layers=None, nodes=None):

        self.node_type = node_type  # cell-centred or vertex-centred
        self.layers = layers
        self.nodes = nodes


class ModelTime(object):
    """
    Class to set the temporal aspects of the model

    within the dictionary "t", the following items are found:
    'start_time': start time of model simulation
    'end_time': end time of model simulation
    'time_step': model time step, if always the same, use 'M' for monthly and 'A' for annual
    'duration': total model simulation time
    'intervals': list of length of time of each time period in the model, stored as datetime object Timedelta
    'steps': number of steps in which model stresses can change
    """

    def __init__(self):
        self.t = {}
        self.t['steady_state'] = True

    def set_temporal_components(self, steady_state=True, start_time=None, end_time=None, time_step=None, date_index=None):
        self.t['steady_state'] = steady_state
        if steady_state == True:
            return

        self.t['start_time'] = start_time
        self.t['end_time'] = end_time
        self.t['time_step'] = time_step
        self.t['duration'] = self.t['end_time'] - self.t['start_time']

        if date_index is None:
            date_index = pd.date_range(start=start_time, end=end_time, freq=time_step)
            self.t['dateindex'] = date_index
        else:
            self.t['dateindex'] = date_index

        intervals = []
        old_time = date_index[0]
        for index, time in enumerate(date_index):
            if index > 0:
                intervals += [time - old_time]
                old_time = time

        self.t['intervals'] = intervals
        self.t['steps'] = len(self.t['intervals'])


class ModelBoundaries(object):
    """
    Class to build boundary conditions with.
    Types can be:
        River
        Recharge
        Wells
        They can be:
        1. At point/s
        2. At layer/s
        3. At domain
        For each type they can be either:
        a. static
        b. dynamic

    Returns a numpy array which for a point will be a location tuple (x, y, z)
    and bc value/s for static, and bc value/s time series for dynamic.
    """

    def __init__(self):
        self.bc = {}
        # self.bc_locale_types = ['point', 'layer', 'domain']
        self.bc_types = ['river', 'wells', 'recharge', 'rainfall',
                         'head', 'drain', 'channel', 'general head']

    def create_model_boundary_condition(self, bc_name, bc_type, bc_static=True, bc_parameter=None):
        # if bc_locale_type not in self.bc_locale_types:
        #    print 'bc_locale_type not recognised, use one of: ', self.bc_locale_types
        #    sys.exit(1)
        # end if
        if bc_type not in self.bc_types:
            print 'bc_type not recognised, use one of: ', self.bc_types
            sys.exit(1)

        self.bc[bc_name] = {}
        self.bc[bc_name]['bc_type'] = bc_type
        self.bc[bc_name]['bc_static'] = bc_static

    def assign_boundary_array(self, bc_name, bc_array):
        if bc_name in self.bc:
            self.bc[bc_name]['bc_array'] = bc_array
        else:
            print 'No boundary condition with name: ', bc_name
            sys.exit(1)


class ModelProperties(object):
    """
    This class is used to set all of the parameters which can then be easily
    accessed for modification for model pertubations
    """

    def __init__(self):
        self.properties = {}
        self.prop_types = ['Kh', 'Kv', 'SS', 'Sy']

    def assign_model_properties(self, prop_type, value):
        if prop_type in self.prop_types:
            self.properties[prop_type] = value
        else:
            print("{} not in {}".format(prop_type, self.prop_types))
            sys.exit('Property type not recognised')


class ModelParameters(object):
    """
    This class is used to set all of the parameters which can then be easily
    accessed for modification for model pertubations
    """

    def __init__(self):
        self.param = {}
        self.param_set = {}

    # 'PARTRANS', 'PARCHGLIM', 'PARVAL1', 'PARLBND', 'PARUBND', 'PARGP', 'SCALE', 'OFFSET'
    def create_model_parameter(self, name, value=None):
        '''
        Function to create new parameter for use in PEST

        '''
        if len(name) > 12:
            print('Warning: PEST has a max char length of parameter names of 12')
            print('         Parameter {0} has length {1}'.format(name, len(name)))
        # end if
        self.param[name] = {}
        self.param[name]['PARVAL1'] = value

    def parameter_options(self, param_name, PARTRANS=None, PARCHGLIM=None,
                          PARLBND=None, PARUBND=None, PARGP=None,
                          SCALE=None, OFFSET=None):
        '''
        Function to assign various paramater properties pertaining to PEST
        '''
        if PARTRANS:
            self.param[param_name]['PARTRANS'] = PARTRANS
        if PARCHGLIM:
            self.param[param_name]['PARCHGLIM'] = PARCHGLIM
        if PARLBND:
            self.param[param_name]['PARLBND'] = PARLBND
        if PARUBND:
            self.param[param_name]['PARUBND'] = PARUBND
        if PARGP:
            self.param[param_name]['PARGP'] = PARGP
        if SCALE:
            self.param[param_name]['SCALE'] = SCALE
        if OFFSET:
            self.param[param_name]['OFFSET'] = PARTRANS

    def create_model_parameter_set(self, name, value=None, num_parameters=1):
        '''
        Function to create a model parameter set to be used with pilot points
        '''
        if len(name) > 9:
            print('Warning: PEST has a max char length of parameter names of 12')
            print('         Parameter {0} has length {1}'.format(name, len(name)))
            print('         Automatic appending of name with number may cause')
            print('         longer than 12 char length par names')
        # end if

        for i in range(num_parameters):
            self.param[name + str(i)] = {}
            self.param[name + str(i)]['PARVAL1'] = value
            if i == 0:
                self.param_set[name] = [name + str(i)]
            else:
                self.param_set[name] += [name + str(i)]

    def parameter_options_set(self, param_set_name, PARTRANS=None, PARCHGLIM=None,
                              PARLBND=None, PARUBND=None, PARGP=None,
                              SCALE=None, OFFSET=None):
        '''
        Function to assign various paramater properties pertaining to PEST
        for each of the parameters within a parameter set
        '''
        for param in self.param_set[param_set_name]:
            self.parameter_options(param,
                                   PARTRANS=PARTRANS, PARCHGLIM=PARCHGLIM,
                                   PARLBND=PARLBND, PARUBND=PARUBND,
                                   PARGP=PARGP, SCALE=SCALE, OFFSET=OFFSET)


class ModelObservations(object):
    """
    This class is used to store all of the obervations relevant to the model
    that are not being used to force the model and for which there are model
    outputs that correspond to the observation, e.g. head
    """

    def __init__(self):
        self.obs_group = {}
        self.obs = {}
        self.obID = 0

    def set_as_observations(self, name, time_series, locations, domain=None,
                            obs_type=None, units=None, weights=None, real=True,
                            by_zone=False):
        """
        Function to set observations from pandas dataframes for times series

        Observations might be: stream stage, stream discharge, stream EC,
        groundwater head, groundwater EC etc.

        Each time series should be of the pandas dataframe format where the
        first column is an identifier for the data, the second column is datetime and
        the next column is the value of interest

        For observation dataframes with multiple identifiers there should be an
        equal number of locations with x, y and z

        """
        self.obs_types = ['head', 'stage', 'discharge', 'concentration']

        self.obs_group[name] = {}
        # check time series meets the standard format of columns = ['name', 'datetime', 'value']
        self.obs_group[name]['time_series'] = time_series
        self.obs_group[name]['time_series']['active'] = True
        self.obs_group[name]['by_zone'] = by_zone
        self.obs_group[name]['time_series']['zone'] = 'null'
        self.obs_group[name]['locations'] = locations
        self.obs_group[name]['domain'] = domain
        self.obs_group[name]['obs_type'] = obs_type
        self.obs_group[name]['units'] = units
        self.obs_group[name]['weights'] = weights
        self.obs_group[name]['real'] = real

    def collate_observations(self):
        for name in self.obs_group.keys():
            self.obs_group[name]['time_series']['obs_map'] = 'null'
            for ob in self.obs_group[name]['time_series'].iterrows():
                if ob[1]['active'] == True:
                    self.obs['ob' + str(self.obID)] = ob[1]['value']
                    self.obs_group[name]['time_series'].set_value(
                        ob[0], 'obs_map', 'ob' + str(self.obID))
                    self.obID += 1

    def check_obs(self):
        obs_nodes = []
        for ob in self.observations.obs.keys():
            obs_nodes += self.observations.obs[ob]


class ModelInitialConditions(object):
    """
    This class is used to store all of the initial conditions for different model
    domains
    """

    def __init__(self):
        self.ic_data = {}

    def set_as_initial_condition(self, name, ic_data):

        self.ic_data[name] = ic_data
        return self.ic_data


class ModelFeature(object):
    """
    This class defines typical features that might be represented in
    a GW model, the definition of which is assigned to particular arrays
    which can then be passed to the appropriate groundwater model.

    """

    def __init__(self, feature_type, feature_name, static=True):

        feature_types = ['Aquifer', 'River', 'Channel', 'Lake',
                         'Wetland', 'Wells', 'Rainfall', 'Evaporation', 'Recharge']
        if feature_type in feature_types:
            self.feature_type = feature_type
        else:
            print 'Feature type ', feature_type, ' not recognised'
            print 'Please use one of: ', feature_types
        self.feature_name = feature_name
        self.static = static


class ModelBuilderType(object):

    """
    This class contains all the types that are allowed for different
    class attributes in GWModelBuilder.
    """

    def __init__(self):
        self.model_types = ['Modflow', 'HGS']
        self.mesh_types = ['structured', 'unstructured']
        self.data_formats = ['ascii', 'binary']
        self.length = ['mm', 'cm', 'm', 'km']
        self.volume = ['ml', 'l', 'kl', 'Ml', 'Gl', 'm3']
        self.time = ['s', 'h', 'd', 'w', 'y']
        self.mass = ['mg', 'g', 'kg']


class ArrayOrdering(object):
    """
    This class describes the array ordering conventions for MODFLOW and HGS
    in structured mesh.

    For the horizontal plane:

    UL  y increasing

        |
        |
        |
        |________
    BL           x increasing

    """

    def __init__(self):
        self.layering_orders = ['TopToBottom', 'BottomToTop']
        self.array_ordering = [
            'UL_RowColumn',  # UL= Upper Left as in MODFLOW
            'BL_RowColumn'   # BL = Bottom Left as in HGS
        ]

    def SetModflowArrays(self):
        self.layer_order = self.layering_orders[0]
        self.array_order = self.array_ordering[0]

    def SetHGSArrays(self):
        self.layer_order = self.layering_orders[1]
        self.array_order = self.array_ordering[1]
