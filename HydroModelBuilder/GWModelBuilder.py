import os
import warnings

import numpy as np
import pandas as pd
from more_itertools import unique_everseen
from scipy import spatial

from ModelInterface.ModelInterface import ModelInterface
from ModelInterface.ModelMesh import MeshGenerator
from ModelProperties import (ArrayOrdering, ModelBoundaries, ModelBuilderType,
                             ModelInitialConditions, ModelObservations,
                             ModelParameters, ModelProperties, ModelTime)
from Utilities import interpolation
from Utilities.PilotPoints import pilotpoints


class GWModelBuilder(object):
    """
    The ModelBuilder class contains a number of useful tools for building
    numerical groundwater models in packages such as MODFLOW, by dealing with
    spatial data using GIS type objects for reading and manipulating spatial
    data before mapping it to a model grid and converting it into an easily readable
    array structure to pass to the groundwater model.

    This class is a first step in the model building process. The object created
    from this class when it is packaged can be used with a separate model interface.
    """

    def __init__(self, name=None, model_type=None, mesh_type=None, units=None,
                 data_folder=None, out_data_folder=None, model_data_folder=None,
                 GISInterface=None, data_format='binary', target_attr=None,
                 **kwargs):
        """
        :param name: str, model name.
        :param model_type: str, Type of model, e.g. MODFLOW (makes use of flopy), HGS (uses pyHGS ... which is not developed yet ;) ).
        :param mesh_type: str, type of mesh to be used in the model, i.e. structured or unstructured.
        :param units:
        :param data_folder: str, Path to get model data from.
        :param out_data_folder: str, Path to store processed model data.
        :param model_data_folder: str, path to model data folder
        :param GISInterface: instance of GISInterface, object to allow interactions with GIS data
        :param data_format: str, format of data files
        :param target_attr: list, target attributes
        """
        # Define the constants for the model data types to use for checking input
        self.types = ModelBuilderType()
        self.types.check_type(model_type, mesh_type, data_format)
        self.name = name
        self.model_type = model_type
        self._mesh_type = mesh_type
        self.ModelInterface = ModelInterface(model_type, self.types, data_format)

        if units is not None:
            if type(units) == list:
                self.set_units(length=units[0], time=units[1], mass=units[2])
            elif type(units) == dict:
                self.set_units(length=units['length'], time=units['time'], mass=units['mass'])
            # End if
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
        # End if

        self.properties = ModelProperties()
        self.boundaries = ModelBoundaries()
        self.parameters = ModelParameters()
        self.observations = ModelObservations()
        self.initial_conditions = ModelInitialConditions()

        # OTHER variables
        self.out_data_folder_grid = None
        self._model_boundary = None
        self.data_boundary = None
        self.boundary_poly_file = None
        self.boundary_data_file = None

        self._model_mesh = None
        self._model_mesh_centroids = None
        self._mesh2centroid2Dindex = None
        self._centroid2mesh2Dindex = None
        self._centroid2mesh3Dindex = None
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

        # Temp object to store river mapping data
        self.river_mapping = {}

        # Object to handle all of the SFR specific data in pandas dataframes
        self.mf_sfr_df = {}

        # Set default target attributes
        if target_attr is None:
            self.target_attr = [
                'name',
                'model_type',
                '_mesh_type',
                '_units',
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
                '_model_boundary',
                'boundary_poly_file',
                'boundary_data_file',
                'model_time',
                '_model_mesh_centroids',
                '_mesh2centroid2Dindex',
                '_model_mesh3D',
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
                'gridWidth',
                'river_mapping',
                'mf_sfr_df'
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
        # End if

        self.MeshGen = MeshGenerator(self.ModelInterface, self.GISInterface)
        self.update_meshgen()
    # End __init__()

    @property
    def model_mesh_centroids(self):
        return self.MeshGen.model_mesh_centroids
    # End model_mesh_centroids()

    @model_mesh_centroids.setter
    def model_mesh_centroids(self, val):
        self._model_mesh_centroids = val
        if hasattr(self, 'MeshGen'):
            self.MeshGen.model_mesh_centroids = val
        # End if
    # End model_mesh_centroids()

    @property
    def mesh_type(self):
        return self.MeshGen.mesh_type
    # End mesh_type()

    @mesh_type.setter
    def mesh_type(self, val):
        self._mesh_type = val
        if hasattr(self, 'MeshGen'):
            self.MeshGen._mesh_type = val
        # End if
    # End mesh_type()

    @property
    def model_mesh(self):
        return self.MeshGen.model_mesh
    # End model_mesh()

    @model_mesh.setter
    def model_mesh(self, val):
        self._model_mesh = val
        if hasattr(self, 'MeshGen'):
            self.MeshGen.model_mesh = val
        # End if
    # End model_mesh()

    @property
    def model_mesh3D(self):
        return self.MeshGen.model_mesh3D
    # End model_mesh3D()

    @model_mesh3D.setter
    def model_mesh3D(self, val):
        self._model_mesh3D = val
        if hasattr(self, 'MeshGen'):
            self.MeshGen.model_mesh = val
        # End if
    # End model_mesh3D()

    @property
    def mesh2centroid2Dindex(self):
        return self.MeshGen.mesh2centroid2Dindex
    # End mesh2centroid2Dindex()

    @mesh2centroid2Dindex.setter
    def mesh2centroid2Dindex(self, val):
        self._mesh2centroid2Dindex = val
        if hasattr(self, 'MeshGen'):
            self.MeshGen.centroid2mesh3Dindex = val
        # End if
    # End mesh2centroid2Dindex()

    @property
    def centroid2mesh2Dindex(self):
        return self.MeshGen.centroid2mesh2Dindex
    # End centroid2mesh2Dindex()

    @centroid2mesh2Dindex.setter
    def centroid2mesh2Dindex(self, val):
        self._centroid2mesh2Dindex = val
        if hasattr(self, 'MeshGen'):
            self.MeshGen.centroid2mesh2Dindex = val
        # End if
    # End centroid2mesh2Dindex()

    @property
    def centroid2mesh3Dindex(self):
        return self.MeshGen.centroid2mesh3Dindex
    # End centroid2mesh2Dindex()

    @centroid2mesh3Dindex.setter
    def centroid2mesh3Dindex(self, val):
        self._centroid2mesh3Dindex = val
        if hasattr(self, 'MeshGen'):
            self.MeshGen.centroid2mesh3Dindex = val
        # End if
    # End centroid2mesh2Dindex()

    @property
    def model_boundary(self):
        return self.MeshGen.model_boundary
    # End model_boundary()

    @model_boundary.setter
    def model_boundary(self, val):
        self._model_boundary = val
        if hasattr(self, 'MeshGen'):
            self.MeshGen.model_boundary = val
        # End if
    # End model_boundary()

    def update_meshgen(self):
        """
        Update MeshGen properties forcibly.
        """
        self.MeshGen.model_boundary = self._model_boundary
        self.MeshGen.model_mesh_centroids = self._model_mesh_centroids
        self.MeshGen.mesh2centroid2Dindex = self._mesh2centroid2Dindex
        self.MeshGen.centroid2mesh2Dindex = self._centroid2mesh2Dindex
        self.MeshGen.centroid2mesh3Dindex = self._centroid2mesh3Dindex
        self.MeshGen.model_mesh3D = self._model_mesh3D
        self.MeshGen.model_mesh = self._model_mesh
        self.MeshGen.mesh_type = self._mesh_type
    # End update_meshgen()

    def update_meshprops(self):
        """
        Forcibly update self with MeshGen properties.
        Can be removed once ModelBuilder and MeshGen are cleanly separated
        """
        self._model_boundary = self.MeshGen.model_boundary
        self._model_mesh_centroids = self.MeshGen.model_mesh_centroids
        self._mesh2centroid2Dindex = self.MeshGen.mesh2centroid2Dindex
        self._centroid2mesh2Dindex = self.MeshGen.centroid2mesh2Dindex
        self._centroid2mesh3Dindex = self.MeshGen.centroid2mesh3Dindex
        self._model_mesh3D = self.MeshGen.model_mesh3D
        self._model_mesh = self.MeshGen.model_mesh
        self._mesh_type = self.MeshGen.mesh_type
    # End update_meshprops()

    def updateGISinterface(self):
        for key, value in self.__dict__.items():
            if type(value) is not object:
                setattr(self.GISInterface, key, value)
            # End if
        # End for
    # End updateGISinterface()

    def set_units(self, length='m', time='s', mass='kg'):
        """
        Sets the units for use inside the model.

        :param length: str, length unit
        :param time: str, time unit
        :param mass: str, mass unit
        """
        self.ModelInterface.set_units(length, time, mass)
        self._units = {"length": length, "time": time, "mass": mass}
    # End set_units()

    def check_for_existing(self, fn):
        """
        Function to determine if input files have previously been processed
        and if so to do nothing unless flagged otherwise. This is done by
        checking the output data path.

        :param fn: str, filename to check for
        """
        self.ModelInterface.check_for_existing(fn)
    # End check_for_existing()

    def save_obj(self, filename, obj):
        self.ModelInterface.save_obj(filename, obj)
    # End save_dataframe()

    def load_obj(self, filename):
        return self.ModelInterface.load_obj(filename)
    # End load_dataframe()

    def save_array(self, filename, array):
        return self.ModelInterface.save_array(filename, array)
    # End save_array()

    def load_array(self, array_file):
        return self.ModelInterface.load_array(array_file)
    # End load_array()

    def save_dataframe(self, filename, df):
        self.ModelInterface.save_dataframe(filename, df)
    # End save_dataframe()

    def load_dataframe(self, filename):
        return self.ModelInterface.load_dataframe(filename)
    # End load_dataframe()

    def flush(self, mode=None):
        self.ModelInterface.flush(mode)
    # End flush()

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
        self.model_boundary[0] = xmin
        self.model_boundary[1] = xmax
        self.model_boundary[2] = ymin
        self.model_boundary[3] = ymax
    # End set_model_boundary_from_corners()

    def set_model_boundary_from_polygon_shapefile(self, shapefile_name, shapefile_path=None):
        if shapefile_path is None:
            shapefile_path = self.data_folder
        # End if

        self.model_boundary, self.boundary_poly_file = self.GISInterface.set_model_boundary_from_polygon_shapefile(
            shapefile_name, shapefile_path)
        return self.model_boundary, self.boundary_poly_file
    # End set_model_boundary_from_polygon_shapefile()

    def set_data_boundary_from_polygon_shapefile(self, shapefile_name=None, shapefile_path=None, buffer_dist=None):
        if shapefile_name is None:
            shapefile_name = self.boundary_poly_file
        if shapefile_path is None:
            shapefile_path = self.out_data_folder
        # End if
        self.data_boundary, self.boundary_data_file = self.GISInterface.set_data_boundary_from_polygon_shapefile(
            shapefile_name, shapefile_path=shapefile_path, buffer_dist=buffer_dist)
        self.updateGISinterface()
        return self.data_boundary, self.boundary_data_file
    # End set_data_boundary_from_polygon_shapefile()

    def build_centroids_array(self, gridHeight):
        """
        1. Builds an array of cell centroids to be used in interpolating from
        other points arrays onto the cell centroids

        2. Creates a dictionary with centroids as key and row col as entries
        If key isn't found then it returns nearest centroid

        :param gridHeight: float, height of grid cells
        """
        return self.MeshGen.build_centroids_array(gridHeight, self.array_ordering)
    # End build_centroids_array()

    def build_centroids_array3D(self):
        """
        1. Builds an array of cell centroids to be used in interpolating from
        other points arrays onto the cell centroids for the 3D mesh

        2. Creates a dictionary with centroids as key and lay row col as entries
        If key isn't found then it returns nearest centroid
        """
        return self.MeshGen.build_centroids_array3D(self.array_ordering)
    # build_centroids_array3D()

    def define_structured_mesh(self, gridHeight, gridWidth):
        """
        Define the structured mesh grid.

        :param gridHeight: float, height of grid cell
        :param gridWidth: float, width of grid cell

        :returns: tuple, grid information
        """
        self.gridHeight = gridHeight
        self.gridWidth = gridWidth
        self.out_data_folder_grid = os.path.join(self.model_data_folder,
                                                 'structured_model_grid_{}m'.format(gridHeight))
        self.updateGISinterface()
        self.model_mesh = self.GISInterface.define_structured_mesh(gridHeight, gridWidth)
        self.model_mesh_centroids = self.build_centroids_array(self.gridHeight)
        return (self.out_data_folder_grid, self.gridHeight, self.gridWidth,
                self.model_mesh_centroids, self.model_mesh)
    # End define_structured_mesh

    def read_rasters(self, files, path=None):
        return self.GISInterface.read_rasters(files, path=path)
    # End read_rasters()

    def map_rasters_to_grid(self, raster_files, raster_path):
        return self.GISInterface.map_rasters_to_grid(raster_files, raster_path)
    # End map_rasters_to_grid()

    def map_raster_to_regular_grid_return_array(self, raster_fname):
        return self.GISInterface.map_raster_to_regular_grid_return_array(raster_fname)
    # End map_raster_to_regular_grid_return_array()

    def map_heads_to_mesh_by_layer(self):
        return self.GISInterface.map_heads_to_mesh_by_layer()
    # End map_heads_to_mesh_by_layer()

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
    # End create_basement_bottom()

    def build_3D_mesh_from_rasters(self, raster_files, raster_path, minimum_thickness, maximum_thickness, force=False):
        """
        :param raster_files: list, raster files to build mesh from
        :param raster_path: str, path to raster files
        :param minimum_thickness: float, minimum thickness
        :param maximum_thickness: float, maximum thickness
        :param force: bool, if True, forces a rebuild of the mesh
        """
        self.MeshGen.build_3D_mesh_from_rasters(self.out_data_folder_grid, self.array_ordering, raster_files,
                                                raster_path, minimum_thickness, maximum_thickness, force)
    # End build_3D_mesh_from_rasters()

    def reclassIsolatedCells(self, passes=1, assimilate=False):
        """
        Function to remove cells that are surrounded by non-active cells in the horizontal plane

        e.g. if cells with positive integer is surrounded in above, below and to each side, then reassign to -1.

        :param passes: int, number of passes to do
        :param assimilate: bool, whether to assimilate or not if surrounded by four of the same.
        """
        self.MeshGen.reclassIsolatedCells(passes, assimilate)
    # End reclassIsolatedCells()

    def read_poly(self, filename, path=None, poly_type='polyline'):
        return self.GISInterface.read_poly(filename, path, poly_type=poly_type)
    # End read_poly()

    def _get_poly_name_and_obj(self, poly):
        """
        Load poly name and object.

        :param poly: str or poly object
        """
        if type(poly) is str:
            poly = self.read_poly(poly)
        # End if
        if os.path.sep in poly.GetDescription():
            poly_name = poly.GetDescription().split(os.path.sep)[-1]
        else:
            poly_name = poly.GetDescription()
        # End if

        return poly_name, poly
    # End _get_poly_name_and_obj()

    def _pointsdist(self, p1, p2, cache={}):
        """
        Calculate distance between points.

        :param p1: array-like, x and y values for point 1
        :param p2: array-like, x and y values for point 2

        :returns: float, distance between points
        """
        tmp = cache.get((p1, p2), None)
        if tmp is not None:
            return tmp
        # End if

        tmp = ((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)**0.5
        cache[(p1, p2)] = tmp

        return tmp
    # End _pointsdist()

    def _hacky_centroid_to_mesh(self, centroid, dist_min):
        """
        Nasty (and slow!) workaround due to minor mismatch in centroids from mesh and separate
        generation in this class. Perhaps better to actually define this array in fishnet when
        define_structured_mesh is called
        """
        return self.MeshGen._hacky_centroid_to_mesh(self.gridHeight, centroid, dist_min, self._pointsdist)
    # End _hacky_centroid_to_mesh()

    ###########################################################################
    # FUNCTIONS USED TO CREATE AND MAP STREAM INFO TO MODEL MESH ... perhaps
    # move these elsewhere but will leave here for now.
    # There is some repitition from above and this needs sorting out!!
    ###########################################################################

    def _dist(self, p1, p2):
        return np.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)
    # End dist()

    def _points_dist_collection(self, points, dist):
        dist_total = 0
        for index in xrange(len(points)):
            if index == 0:
                pass
            else:
                dist_total += dist(points[index], points[index - 1])
        return dist_total
    # End points_dist_collection()
    
    def create_river_dataframe(self, name, poly_file, surface_raster_file,
                               plotting=False, avoid_collocation=False, verbose=False):

        poly_file = self.GISInterface.read_poly(poly_file)
        points_list, points_dict = self.GISInterface.polyline_explore(poly_file)
        points_dict = self._get_start_and_end_points_from_line_features(points_dict)
        point_merge = self._points_merge(points_dict)
        if verbose:
            print("Find closest model mesh cell for line points")
        closest = self.do_kdtree(np.array(self.centroid2mesh2Dindex.keys()), point_merge)
        point2mesh_map, point2mesh_map2 = self._points2mesh_map(point_merge, closest)
        amalg_riv_points, amalg_riv_points_collection = self._amalgamate_points(point2mesh_map2, point_merge)
        # Do test and report any cell jumps more than 1 up, down, left, right
        if verbose:
            print("Test order the cells")
        self._naive_cell_ordering_test(amalg_riv_points)

        dist = self._dist

        if verbose:
            print("Get length of reaches")
        lengths = []
        for i in xrange(len(amalg_riv_points_collection)):
            if i == 0:
                lengths += [self._points_dist_collection(amalg_riv_points_collection[i], dist)]
                continue
            # End if
            lengths += [self._points_dist_collection([amalg_riv_points_collection[i - 1][-1]] +
                                               amalg_riv_points_collection[i], dist)]
        # End for

        if verbose:
            print("Process the surface raster for elevations")
            
        surf_raster_fname = os.path.join(self.out_data_folder_grid, 'surf_raster_processed_{}.pkl'.format(name))
        if os.path.exists(surf_raster_fname):
            print(" --- Using previously processed surface raster data --- ")
            sr_array, point2surf_raster_map = self.ModelInterface.load_obj(surf_raster_fname)
        else:
            surf_raster_info = self.get_raster_info(surface_raster_file)
            srm = surf_raster_info['metadata']
            if verbose:
                print("Calculate centroids")
            surf_centroids = self._create_centroids(srm['pixel_x'], srm['pixel_y'], (srm['ulx'],
                                                                                     srm['ulx'] + srm['cols'] *
                                                                                     srm['pixel_x'],
                                                                                     srm['uly'] - srm['rows'] *
                                                                                     srm['pixel_y'],
                                                                                     srm['uly']))
            
            surf_raster_points = np.array(surf_centroids.keys())
            if verbose:
                print("Using kd_tree to find surface raster cells closest to river points")
            closest_srp = self.do_kdtree(surf_raster_points, point_merge)
            point2surf_raster_map = []

            if verbose:
                print("creating the point2surf_raster_map")
            for index, point in enumerate(point_merge):
                point2surf_raster_map += [surf_centroids[tuple(surf_raster_points[closest_srp[index]])]]
            sr_array = surf_raster_info['array']
            self.ModelInterface.save_obj([sr_array, point2surf_raster_map], surf_raster_fname[:-4])
        # End if

        if verbose:
            print("Assign elevations based on surface raster map")

        riv_elevations = [sr_array[x[0]][x[1]] for x in point2surf_raster_map]
        riv_reach = []
        reach = 0

        for index, point in enumerate(point2mesh_map2):
            # If it is the last point
            if index == len(point2mesh_map2) - 1:
                riv_reach += [reach]
                continue
            if index == 0:
                riv_reach += [reach]
                continue
            if point == point2mesh_map2[index + 1]:
                riv_reach += [reach]
            elif point != point2mesh_map2[index + 1]:
                reach += 1
                riv_reach += [reach]
            # End if
        # End for

        if verbose:
            print("Creating river dataframe")

        river_df = pd.DataFrame({'reach': riv_reach, 'strtop': riv_elevations})
        # In case of hitting no data values replace these with nan
#        river_df.loc[river_df['strtop'] < 0, 'strtop'] = np.nan
        # Then to fill the Nan values interpolate or drop??
#        river_df.dropna(inplace=True)

        # Group the dataframe stream points by reach
        river_seg = river_df.groupby('reach').min()
        river_seg.loc[:, 'rchlen'] = lengths
        river_seg.loc[:, 'amalg_riv_points'] = amalg_riv_points

        # Check that order of elevations if from upstream to downstream ...
        elev = river_seg['strtop'].tolist()
        # Reorder the dataframe if this is the case
        if elev[0] < elev[-1]:
            river_seg = river_seg.reindex(index=river_seg.index[::-1])
        # End if

        # Now get the cumulative length along the stream
        river_seg['Cumulative Length'] = river_seg['rchlen'].cumsum()

        river_reach_elev = river_seg['strtop'].tolist()
        river_reach_elev_adjusted = []
        for index, elev_reach in enumerate(river_reach_elev):
            if index == 0:
                river_reach_elev_adjusted += [elev_reach]
                continue
            # End if
            if index == len(river_reach_elev) - 1:
                river_reach_elev_adjusted += [self._adjust_elevations(
                    river_reach_elev_adjusted[index - 1], elev_reach, None)]
                continue
            # End if

            river_reach_elev_adjusted += [self._adjust_elevations(
                river_reach_elev_adjusted[-1], elev_reach, river_reach_elev[index + 1])]
        # End for
        river_seg.loc[:, 'strtop_raw'] = river_seg['strtop']
        river_seg.loc[:, 'strtop'] = river_reach_elev_adjusted
        if plotting:
            ax = river_seg.plot(x='Cumulative Length', y='strtop_raw', alpha=0.3)
            river_seg.plot(x='Cumulative Length', y='strtop', ax=ax)
        # End if

        slopes = []
        len_rr_elev_adjusted = len(river_reach_elev_adjusted) - 1
        for index, riv_elev in enumerate(river_reach_elev_adjusted):
            if index == len_rr_elev_adjusted:
                slopes += [slopes[-1]]
                continue
            # End if
            slopes += [(riv_elev - river_reach_elev_adjusted[index + 1]) / lengths[index]]
        # End for

        river_seg.loc[:, 'slope'] = slopes

        amalg_riv_points_naive_layer = [[0] + x for x in river_seg['amalg_riv_points'].tolist()]
        river_seg.loc[:, 'k'] = [x[0] for x in amalg_riv_points_naive_layer]
        river_seg.loc[:, 'i'] = [x[1] for x in amalg_riv_points_naive_layer]
        river_seg.loc[:, 'j'] = [x[2] for x in amalg_riv_points_naive_layer]
        river_seg.loc[:, 'amalg_riv_points_collection'] = [amalg_riv_points_collection[x]
                                                    for x in xrange(len(amalg_riv_points_collection.keys()))]

        self.river_mapping[name] = river_seg

    def get_closest_riv_segments(self, name, points):
        '''
        Fuction to find the closest river segment defined in self.river_mapping
        for the points passed to the function.

        This function requires that self.river_mapping exists, which occurs
        after running "self.create_river_dataframe"

        :param name: str, Name given to stream or river when create_river_dataframe was run
        :param points: list, points for which to find the closest to named river

        :returns: list, closest river segment
        '''
        try:
            self.river_mapping[name]
        except NameError:
            print("This function can only be run after running 'create_river_dataframe'")
        # End except

        river_points = []
        amalg_riv_points = self.river_mapping[name]['amalg_riv_points_collection'].tolist()
        riv_len = self.river_mapping[name].shape[0]
        for i in xrange(riv_len):
            river_points += amalg_riv_points[i]
        # End for
        closest = self.do_kdtree(np.array(river_points), np.array(points))

        river_seg_close = []
        for index in closest:
            for i in xrange(riv_len):
                if river_points[index] in amalg_riv_points[i]:
                    river_seg_close += [i + 1]  # Note that seg is 1-based indexing so we add 1 here
                    continue
                # End if
            # End for
        # End for
        return river_seg_close
    # End get_closest_riv_segments()

    def _adjust_elevations(self, us, active, ds, adjust=0.1, verbose=False):
        '''
        Function to force the river to be decreasing in elevation as you move downstream

        :param us: float, upstream elevation
        :param active: float, elevation being analysed
        :param ds: float, downstream elevation
        :param adjust: float, amount to adjust by
        :param verbose: bool, print out debug information
        '''
        if us is None:
            if verbose:
                print("No upstream to worry about, move on")
            return active

        if us > active and active > ds:
            if verbose:
                print("The active elevation fits between upstream and downstream, move on")
            return active

        if us < active:
            if ds is not None:
                if us > ds:
                    if verbose:
                        print("Active greater than upstream and downstream, interpolate between us and ds")
                    return (us + ds) / 2.
                if us < ds and verbose:
                    print("Case 4")
            else:
                if verbose:
                    print("Case 5")
            # End if
        # End if

        if us == active and verbose:
            print("Case 6")

        if verbose:
            print("Case 7: {} {} {}".format(us, active, ds))
        return us - adjust
    # End _adjust_elevations()

    def _create_centroids(self, x_pixel, y_pixel, bounds):
        '''
        Function to create centroids of raster cells given the bounds of the
        raster and the width and height of the cells

        :param x_pixel: float,
        :param y_pixel: float,
        :param bounds: tuple,

        :returns: tuple, (dict) centroid2mesh, (dict) mesh2centroid
        '''
        xmin, xmax, ymin, ymax = bounds
        print bounds, x_pixel, y_pixel
        cols = int((xmax - xmin) / x_pixel)
        rows = int((ymax - ymin) / y_pixel)
        x = np.linspace(xmin + x_pixel / 2.0, xmax - x_pixel / 2.0, cols)
        y = np.linspace(ymax - y_pixel / 2.0, ymin + y_pixel / 2.0, rows)
        
        centroid2mesh2Dindex = {}
        for row in xrange(rows):
            for col in xrange(cols):
                centroid2mesh2Dindex[(x[col], y[row])] = [row, col]
            # End for
        # End for
        print "finished centroid_building"
        return centroid2mesh2Dindex
    # End _create_centroids()

    def get_raster_info(self, raster_file):
        '''
        Function to get raster data and array into python objects

        :param raster_file: str,
        '''
        return self.GISInterface.get_raster_info(raster_file)
    # End get_raster_info()

    def _get_lengths_for_polyline_in_cells(self, point_merge, point2mesh_map2):
        """
        :param point_merge:
        :param point2mesh_map2:
        """
        lengths = []
        current = 0.0
        carryover = 0.0

        dist = self._dist
        for index, points in enumerate(point2mesh_map2):
            index_next = index + 1
            # Stop at the last iteration
            if index == len(point2mesh_map2) - 1:
                lengths += [current + carryover]
                continue
            # If points are in same cell, calculage length and add to current
            if point2mesh_map2[index] == point2mesh_map2[index_next]:
                current += dist(point_merge[index], point_merge[index_next])
            elif point2mesh_map2[index] != point2mesh_map2[index_next]:
                lengths += [current + carryover]
                current = 0.
                carryover = dist(point_merge[index], point_merge[index_next])
            # End if
        # End for
    # End _get_lengths_for_polyline_in_cells()

    def _naive_cell_ordering_test(self, cell_list):
        '''
        Function to test the ordering of cells and see that they are sequentially
        neighbours and not jumping multiple cells.

        This is a test for river cell ordering which is specific for MODFLOW
        SFR and STR packages.

        :param cell_list: list, of cells to order
        '''
        cl = cell_list
        for index in xrange(len(cell_list)):
            if index == 0:
                pass
            else:
                if np.sqrt((cl[index][0] - cl[index - 1][0]) ** 2 + (cl[index][1] - cl[index - 1][1]) ** 2) > 1:
                    print("Warning, cell jump from {} to {}".format(cl[index - 1], cl[index]))
                # End if
            # End if
        # End for
    # _naive_cell_ordering_test()

    def _amalgamate_points(self, point2mesh_map2, point_merge):
        """
        Amalgamate points.

        :param point2mesh_map2:
        :param point_merge:

        :returns: tuple, (list) river points, (dict) river points collection
        """
        # TO CHECK
        amalg_riv_points = []
        amalg_riv_points_collection = {}
        counter = 0
        for index, point in enumerate(point2mesh_map2):
            amalg_riv_points_collection[counter] = amalg_riv_points_collection.get(counter, [])
            if index == 0:
                amalg_riv_points += [point]
            elif point != amalg_riv_points[-1]:
                # if point not in amalg_riv_points:  # Apply this test to avoid passing through same cell twice!
                amalg_riv_points += [point]
                counter += 1
                amalg_riv_points_collection[counter] = amalg_riv_points_collection.get(counter, [])
            # End if

            amalg_riv_points_collection[counter] += [point_merge[index]]
        # End for

        return amalg_riv_points, amalg_riv_points_collection
    # End _amalgamate_points()

    def _points2mesh_map(self, points, closest):
        """
        :param points:
        :param closest:

        :returns: tuple, (dict) point2mesh mapping, (list) point2mesh in list form
        """
        points_array = np.array(points)
        centroids = self.centroid2mesh2Dindex
        model_mesh_points = np.array(centroids.keys())

        point2mesh_map = {}
        for index, point in enumerate(points_array):
            point2mesh_map[tuple(list(point))] = centroids[tuple(model_mesh_points[closest[index]])]
        # End for

        warnings.warn("This method will only return a dictionary in future. List `point2mesh_map2` will be removed",
                      DeprecationWarning)
        point2mesh_map2 = []
        for index, point in enumerate(points_array):
            point2mesh_map2 += [centroids[tuple(model_mesh_points[closest[index]])]]
        # End for

        return point2mesh_map, point2mesh_map2
    # End _points2mesh_map()

    def _get_start_and_end_points_from_line_features(self, points_dict):
        """
        :param points_dict:
        """
        for key in points_dict:
            points_dict[key] += [{'start': points_dict[key][0], 'end':points_dict[key][-1]}]
        return points_dict
    # End _get_start_and_end_points_from_line_features()

    def _points_merge(self, points_dict):
        '''
        Merge points from multiple feature into one continuous line

        :param points_dict: dict, of feature points

        :return: list, of unique points.
        '''
        # TO CHECK
        point_merge = points_dict[0][0:-1]
        points_after_first = points_dict.keys()[1:]
        for _ in points_after_first:
            for key in points_after_first:
                if points_dict[key][-1]['start'] == point_merge[-1]:
                    point_merge += points_dict[key][0:-1]
                    continue
                # End if
                if points_dict[key][-1]['end'] == point_merge[0]:
                    point_merge = points_dict[key][0:-1] + point_merge
                    continue
                # End if
            # End for
        # End for
        point_merge = list(unique_everseen(point_merge))

        return point_merge
    # End _points_merge()

    def save_MODFLOW_SFR_dataframes(self, name, reach_df, seg_df):
        self.mf_sfr_df[name] = {'reach_df': reach_df, 'seg_df': seg_df}
    # End save_MODFLOW_SFR_dataframes()

    ###########################################################################
    ###########################################################################
    ###########################################################################

    def map_polyline_to_grid(self, polyline_obj):
        poly_name, polyline_obj = self._get_poly_name_and_obj(polyline_obj)
        if os.path.exists(os.path.join(self.out_data_folder_grid, poly_name + '_mapped.pkl')):
            print "Using previously mapped polyline to grid object"
            self.polyline_mapped[poly_name] = self.ModelInterface.load_obj(os.path.join(self.out_data_folder_grid,
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
                    grid_loc = self._hacky_centroid_to_mesh(centroid, dist_min=1E6)

                temp += [[grid_loc, item[0]]]
            # End for
            self.polyline_mapped[poly_name] = temp
            self.ModelInterface.save_obj(temp, os.path.join(self.out_data_folder_grid, poly_name + '_mapped'))
        # End if
        self.gridded_data_register += [poly_name]
    # End map_polyline_to_grid()

    def map_polygon_to_grid(self, polygon_obj, feature_name=None):

        poly_name, polygon_obj = self._get_poly_name_and_obj(polygon_obj)
        if os.path.exists(os.path.join(self.out_data_folder_grid, poly_name + '_mapped.pkl')):
            print "Using previously mapped polygons to grid object"
            self.polygons_mapped[poly_name] = self.ModelInterface.load_obj(os.path.join(self.out_data_folder_grid,
                                                                                        poly_name + '_mapped.pkl'))
        else:
            self.polygons_mapped[poly_name] = \
                self.GISInterface.map_polygon_to_grid(polygon_obj,
                                                      out_fname=os.path.join(
                                                          self.out_data_folder_grid, poly_name + '_mapped'),
                                                      pixel_size=self.model_mesh_centroids[0][0][
                                                          1] - self.model_mesh_centroids[0][0][0],
                                                      bounds=self.model_mesh.GetLayer().GetExtent(),
                                                      feature_name=feature_name)

            self.ModelInterface.save_obj(self.polygons_mapped[poly_name],
                                         os.path.join(self.out_data_folder_grid,
                                                      poly_name + '_mapped'))
        # End if

        self.gridded_data_register += [poly_name]
    # End map_polygon_to_grid()

    def read_points_data(self, filename, path=None):
        return self.GISInterface.read_points_data(filename, path)
    # End read_points_data()

    def map_points_to_grid(self, points_obj, feature_id=None, use_kdtree=False):
        point_name, points_obj = self._get_poly_name_and_obj(points_obj)
        if os.path.exists(os.path.join(self.out_data_folder_grid, point_name + '_mapped.pkl')):
            print "Using previously mapped points to grid object"
            self.points_mapped[point_name] = self.ModelInterface.load_obj(os.path.join(self.out_data_folder_grid,
                                                                                       point_name + '_mapped.pkl'))
        else:
            temp = []
            if use_kdtree:
                # Need to get points from points_obj
                points = []
                points_label = []
                
                # Get feature id's as well from 
                identifier = ""
                temp = self.map_points_to_2Dmesh(points, identifier)
            else:
                self.points_mapped[point_name] = self.GISInterface.map_points_to_grid(
                    points_obj, feature_id=feature_id)
                # Go from centroids to ij indexing
                for item in self.points_mapped[point_name]:
                    centroid = (float(item[1][0]), float(item[1][1]))
                    # replace centroid with row col
                    try:
                        grid_loc = self.centroid2mesh2Dindex[centroid]
                    except:
                        grid_loc = self._hacky_centroid_to_mesh(centroid, dist_min=1E6)

                    temp += [[grid_loc, item[0]]]

            self.points_mapped[point_name] = temp
            self.ModelInterface.save_obj(temp, os.path.join(self.out_data_folder_grid, point_name + '_mapped'))

        self.gridded_data_register += [point_name]
    # End map_points_to_grid()

    def do_kdtree(self, model_mesh_points, points):
        mytree = spatial.cKDTree(model_mesh_points)
        dist, indexes = mytree.query(points)
        return indexes
    # End do_kdtree()

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
        # End if

        closest = self.do_kdtree(model_mesh_points, points)
        point2mesh_map = {}
        for index, point in enumerate(points):
            if identifier[0]:
                point2mesh_map[identifier[index]] = self.centroid2mesh2Dindex[
                    tuple(model_mesh_points[closest[index]])]
            else:
                point2mesh_map[point] = self.centroid2mesh2Dindex[
                    tuple(model_mesh_points[closest[index]])]
            # End if
        # End for

        return point2mesh_map
    # End map_points_to_2Dmesh()

    def find_closest_points_between_two_lists(self, list1, list2):
        np_arr1 = np.array(list1)
        np_arr2 = np.array(list2)
        closest = self.do_kdtree(np_arr1, np_arr2)

        return closest
    # End find_closest_points_between_two_lists()

    def map_points_to_3Dmesh(self, points, identifier=None):
        '''
        Function to map points to the 3D mesh

        Returns: A dict with keys of all points (or identifiers) with each
        corresponding entry containing the i, j and k reference of the nearest
        cell center to the given point
        '''
        return self.MeshGen.map_points_to_3Dmesh(self.do_kdtree, points, identifier)
    # End map_points_to_3Dmesh()

    def map_obs_loc2mesh3D(self, method='nearest', ignore=[-1]):
        """
        This is a function to map the obs locations to the nearest node in the
        mesh
        """
        self.MeshGen.map_obs_loc2mesh3D(self.observations, self.do_kdtree, method, ignore)
    # End map_obs_loc2mesh3D()

    def getXYpairs(self, points_obj, feature_id=None):
        """
        TODO: docs
        """
        return self.GISInterface.getXYpairs(points_obj, feature_id=feature_id)
    # End getXYpairs()

    def points_shapefile_obj2dataframe(self, points_obj, feature_id=None):
        """
        Converts GIS object to dataframe, which is designed for point style objects only
        """
        points_dict = self.GISInterface.getXYpairs(points_obj, feature_id=feature_id)
        df_points = pd.DataFrame.from_dict(points_dict, orient='index')
        df_points.columns = ['Easting', 'Northing']
        df_points[feature_id] = df_points.index

        return df_points
    # End points_shapefile_obj2dataframe()

    def point_values_from_raster(self, points, raster_obj):
        """
        Get values from raster at given points defined by list or by points object
        points must be a list made up of point lists or points shapefile object
        """
        return self.GISInterface.point_values_from_raster(points, raster_obj)
    # End point_values_from_raster()

    def map_points_to_raster_layers(self, points, depths, rasters):
        """
        TODO: docs
        """
        # create boolean array
        len_list = len(points)
        layers = len(rasters) / 2  # Intentional integer division
        points_layer = np.full((layers, len_list), False, dtype=bool)

        return self.GISInterface.map_points_to_raster_layers(points, depths, rasters, points_layer)
    # End map_points_to_raster_layers()

    def interpolate_points2mesh(self, points_obj, values_dataframe, feature_id=None,
                                method='nearest', use='griddata',
                                function='multiquadric', epsilon=2):
        """
        Interpolate points to mesh.

        :param points_obj: array-like, points to interpolate
        :param values_dataframe: DataFrame, of values
        :param feature_id: int, feature id defaulting to `None`
        :param method: str, interpolation method. Defaults to `nearest`
        :param use: str, Defaults to `griddata`
        :param function: str, Defaults to `multiquadric`
        :param epsilon: int, Defaults to `2`
        """
        if isinstance(points_obj, list) or isinstance(points_obj, np.ndarray):
            points = points_obj
            values = values_dataframe
        else:
            points_dict = self.getXYpairs(points_obj, feature_id=feature_id)
            points = []
            values = []
            for key in points_dict:
                try:
                    values += [float(values_dataframe[key])]
                    points += [points_dict[key]]
                except:
                    pass

        return interpolation.Interpolator(self.mesh_type, np.array(points), np.array(values),
                                          self.model_mesh_centroids, method=method, use=use,
                                          function=function, epsilon=epsilon)
    # End interpolate_points2mesh()

    @staticmethod
    def _findInterval(row, times):
        """
        TODO: docs
        """
        key_time = row['datetime']
        lower_time = times[0]
        for period, time in enumerate(times):
            if period > 0:
                if lower_time <= key_time < time:
                    return period - 1
                # End if
            # End if
            lower_time = time
        # End for

        return np.nan
    # End _findInterval()

    def map_obs2model_times(self):
        """
        This is a function to map the obs at different times within the bounding interval in the
        model times intervals
        """
        obs_group = self.observations.obs_group
        if self.model_time.t['steady_state']:
            for key in obs_group:
                obs_group[key]['time_series']['interval'] = 0
        else:
            find_interval = self._findInterval
            dateindex = self.model_time.t['dateindex']
            for key in obs_group:
                obs_group_key = obs_group[key]
                obs_group_key['time_series'].loc[:, 'interval'] = obs_group_key['time_series'].apply(
                    lambda row: find_interval(row, dateindex), axis=1)
                # remove np.nan values from the obs as they are not relevant
                if obs_group_key['real']:
                    obs_group_ts = obs_group_key['time_series']
                    obs_group_ts = obs_group_ts[pd.notnull(obs_group_ts['interval'])]
                    obs_group_ts = obs_group_ts[obs_group_ts['value'] != 0.]
                # End if
            # End for
        # End if
    # End map_obs2model_times()

    def updateModelParameters(self, fname, verbose=True):
        """
        Update model parameters based on values given in a file.

        :param fname: str, filename to open
        :param verbose: bool, whether or not to print out debug info
        """
        with open(fname, 'r') as f:
            text = f.readlines()
            # Remove header
            text = text[1:]
            # Read in parameters and replace values in parameters class for param
            updated = {}
            for key in self.parameters.param:
                updated[key] = False

            for line in text:
                param_name, value = line.strip('\n').split('\t')
                value = value.lstrip()
                param_name = param_name.strip()
                if param_name in self.parameters.param:
                    self.parameters.param[param_name]['PARVAL1'] = float(value)
                    updated[param_name] = True
                else:
                    if verbose:
                        print 'Parameter not defined in model: ', param_name

        if verbose:
            were_updated = [key for key in updated if updated[key] == True]
            if len(were_updated) > 0:
                print 'Parameters updated for : ', were_updated
            not_updated = [key for key in updated if updated[key] == False]
            if len(not_updated) > 0:
                print 'Parameters unchanged for : ', not_updated
    # End updateModelParameters()

    def generate_update_report(self):
        """"
        Generate report on boundaries and properties that have been updated

        For use in model running to report on which boundaries and properties
        have been changed in the run scripts. Requires use of update commands
        for both properties and boundaries to work, i.e.,
            boundaries.update_boundary_array
            update_model_properties
        """
        self.boundaries.generate_update_report()
        self.properties.generate_update_report()
    # End generate_update_report()

    def create_pilot_points(self, name, linear=False):
        """
        :param name: str, identifier of target pilot point
        :param linear: bool, linear or non-linear pilot points.
        """
        if not linear:
            self.pilot_points[name] = pilotpoints.PilotPoints(
                output_directory=self.out_data_folder_grid, additional_name=name)
        else:
            self.pilot_points[name] = pilotpoints.PilotPointsLinear()
        # End if
    # End create_pilot_points()

    def save_pilot_points(self):
        """
        Save pilot point data in pickle file.
        """
        self.ModelInterface.save_obj(self.pilot_points, os.path.join(self.out_data_folder_grid, 'pilot_points'))
    # End save_pilot_points()

    def load_pilot_points(self, fname):
        """
        Load pilot point data from pickle file.

        :param fname: str, filename to load data from.
        """
        pp = self.ModelInterface.load_obj(fname)
        for key in pp:
            self.pilot_points[key] = pp[key]
        # End for
    # End load_pilot_points()

    def add2register(self, addition):
        """
        Add given value to model register.

        :param addition: numeric, value to add.
        """
        warnings.warn("Call to (possibly) deprecated method `add2register`", DeprecationWarning)
        self.model_register += addition
    # End add2register()

    def writeRegister2file(self):
        """
        Write model register to file.
        """
        with open(os.path.join(self.out_data_folder, 'model_register.dat')) as f:
            for item in self.model_register:
                f.write(item)
            # End for
        # End with
    # End writeRegister2file()

    def points2shapefile(self, points_array, shapefile_name):
        """
        Write given points array to shapefile.

        :param points_array:
        :param shapefile_name:
        """
        warnings.warn("Call to (possibly) deprecated method `points2shapefile`", DeprecationWarning)
        self.GISInterface.points2shapefile(points_array, shapefile_name)
    # End points2shapefile()

    def mesh3DToVtk(self, val_array, val_name, out_path, vtk_out):
        '''
        Function to write the mesh array.

        :param val_array:
        :param val_name:
        :param out_path:
        :param vtk_out:
        '''
        from HydroModelBuilder.GISInterface.GDALInterface import array2Vtk
        warnings.warn("Call to (possibly) deprecated method `mesh3DToVtk`", DeprecationWarning)
        nrow, ncol = self.model_mesh3D[1][0].shape
        delc, delr = self.gridWidth, self.gridHeight
        x0, y0 = self.model_boundary[0], self.model_boundary[3]
        grid_info = [ncol, nrow, delc, delr, x0, y0]
        mesh, zone_matrix = self.mesh_array[0], self.mesh_array[1]

        array2Vtk.build_vtk_from_array(grid_info, np.fliplr(mesh), ["z_elev"],
                                       [np.fliplr(mesh)], ["zone", val_name],
                                       [np.fliplr(zone_matrix), np.fliplr(val_array)],
                                       out_path, vtk_out)
    # End mesh3DToVtk()

    def package_data(self):
        """
        Option to save all important attributes of GWModelBuilder class to
        allow quick loading of data that may have required transforms and
        processing in its orignial state.

        This will include GIS type objects
        """
        pass
    # End package_data()

    def package_model(self):
        """
        Option to save all important attributes of model class to allow quick
        loading and manipulation of model without separate data and grid generation
        building scripts

        This will not include any GIS type objects
        """
        self.update_meshprops()
        self.ModelInterface.package_model(self.out_data_folder_grid, self.name + '_packaged',
                                          self.__dict__, self.target_attr)
    # End package_model()

    ###########################################################################

    def __exit__(self):
        self.writeRegister2file()
        self.save_mapped_dictionaries()
    # End __exit__()
# End GWModelBuilder()
