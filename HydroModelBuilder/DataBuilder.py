# built-ins
import sys
import os
import shutil
import pickle

# additions
import numpy as np
import pandas as pd

class DataBuilder(object):
    """
    The DataBuilder class contains a number of useful tools for collecting
    temporal and spatial datasets and bringing them into one place for then
    building numerical groundwater models in packages such as MODFLOW.
    """

    def __init__(self, name=None, units=None, data_folder=None, out_data_folder=None,
                 GISInterface=None, data_format='binary', target_attr=None, **kwargs):

        """
        :param name: Dataset name.
        :param data_folder: Path to get model data from.
        :param out_data_folder: Path to store processed model data.


        """
        # Define the constants for the model data types to use for checking input
        self.types = ModelBuilderType()

        try:
            # -- tests to alert user to incorrect inputs ...
            assert data_format in self.types.data_formats, "Data format must be of type: {}".format(
                self.types.data_formats)

            # if data_folder != None:
            #    assert os.path.isdir(data_folder) == True, "{} is an invalid path".format(data_folder)
            # End if

            # assert os.path.isdir(out_data_folder) == True, "{} is an invalid path".format(out_data_folder)

        except AssertionError as e:
            import traceback
            _, _, tb = sys.exc_info()
            # traceback.print_tb(tb) # Fixed format
            tb_info = traceback.extract_tb(tb)
            filename, line, func, text = tb_info[-1]
            sys.exit("An error occured in {} on line {} with the message '{}'".format(filename, line, e))
        # End try

        self.name = name

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

        # OTHER variables
        self.out_data_folder_grid = None
        self.model_boundary = None
        self.data_boundary = None
        self.boundary_poly_file = None
        self.boundary_data_file = None

        # Create registers to store details of processed and imported data for
        # quick future loading
        self.model_register = []
        self.base_data_register = []
        self.gridded_data_register = []

        # Set default target attributes
        if target_attr is None:
            self.target_attr = [
                'name',
                'units',
                'data_folder',
                'out_data_folder',
                'data_format',
                'array_ordering',
                'boundaries',
                'properties',
                'parameters',
                'initial_conditions',
                'out_data_folder_grid',
                'model_boundary',
                'boundary_poly_file',
                'boundary_data_file',
                'base_data_register',
            ]
        else:
            self.target_attr = target_attr
        # End if

        # Set all other kwargs as class attributes
        for key, value in kwargs.items():
            setattr(self, key, value)
        # End For

        # Set all key word arguments as attributes
        if self.GISInterface == None:
            print 'No GIS interface defined'
        else:
            for key, value in self.__dict__.items():
                if type(value) is not object:
                    setattr(self.GISInterface, key, value)
                    #setattr(self.ModelInterface, key, value)
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

    ###########################################################################
    # START GENERIC FUNCTIONS WHICH SHOULD BE IMPORTED    
    ###########################################################################

    def check_for_existing(self, f):
        """
        Function to determine if input files have previously been processed
        and if so to do nothing unless flagged otherwise. This is done by
        checking the output data path.
        """
        filename_suffixes = ['_model', '_grid']

        for suffix in filename_suffixes:
            if os.path.isfile(self.out_data_folder + f[:-4] +
                              suffix + f[-4:]):
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
        if array_file[-3:] == 'txt':
            return np.loadtxt(array_file)
        elif array_file[-3:] == 'npy' or array_file[-3:] == 'npz':
            return np.load(array_file)
        else:
            print 'File type not recognised as "txt", "npy" or "npz" \n'
            sys.exit(1)
        # end if

    def save_dataframe(self, filename, df):
        df.to_hdf(filename + '.h5', 'table')

    def load_dataframe(self, filename):
        if filename[-3:] == '.h5':
            return pd.read_hdf(filename, 'table')
        else:
            print 'File type not recognised as "h5"'
            sys.exit(1)
        # end if

    def save_obj(self, obj, filename):
        with open(filename + '.pkl', 'wb') as f:
            pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

    def load_obj(self, filename):
        if filename[-4:] == '.pkl':
            with open(filename, 'rb') as f:
                return pickle.load(f)
        else:
            print 'File type not recognised as "pkl"'
            sys.exit(1)
        # end if

    def flush(self):
        folder = self.out_data_folder
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

    ###########################################################################
    # END GENERIC FUNCTIONS WHICH SHOULD BE IMPORTED    
    ###########################################################################

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

    def set_model_boundary_from_polygon_shapefile(self, shapefile_name, shapefile_path=None):
        """
        Function to set the model boundary from a polygonal shapefile
        Assumes that there is only one polygon in the shapefile!
        
        :param shapefile_name: Name of the shapefile including extension (String) [required]
        :param shapefile_path: Path for the shapefile is not current working directory (String) [optional]
        """
        self.model_boundary, self.boundary_poly_file = self.GISInterface.set_model_boundary_from_polygon_shapefile(
            shapefile_name, shapefile_path)
        return self.model_boundary, self.boundary_poly_file

    def set_data_boundary_from_polygon_shapefile(self, shapefile_name, shapefile_path=None, buffer_dist=None):
        """
        Function to set the data boundary from a polygonal shapefile which also
        allows use of a buffer distance around the polygon so that if the same 
        shapefile that is used to set the model is used, the data boundary can
        encompass a wider reaching area.

        The data boundary allows for filtering of large spatial data sets to only
        include data within the data boundary. This is useful whereby the data 
        within a specified model area as well as within some distance outside of
        the model area can be used to inform the model.
        
        Assumes that there is only one polygon in the shapefile!
        
        :param shapefile_name: Name of the shapefile including extension (String) [required]
        :param shapefile_path: Path for the shapefile is not current working directory (String) [optional]
        :param buffer_dist: Buffer distance around polygon to expand coverage (Float) [optional]
        """

        self.data_boundary, self.boundary_data_file = self.GISInterface.set_data_boundary_from_polygon_shapefile(
            shapefile_name, shapefile_path=shapefile_path, buffer_dist=buffer_dist)
        self.updateGISinterface()
        return self.data_boundary, self.boundary_data_file

    def read_rasters(self, files, path=None):
        """
        Function to read in rasters for use elsewhere
        """
        return self.GISInterface.read_rasters(files, path=path)
        
        
    def create_basement_bottom(self, hu_raster_path, surface_raster_file, basement_top_raster_file, basement_bot_raster_file, output_path, raster_driver='GTiff'):
        """
        Utility to build a bottom basement array where it doesn't exist based on top of bedrock, surface elevation and a thickness function

        writes a raster file for the bottom basement array

        ** This perhaps would better sit in a separate utilities folder ...
        """
        return self.GISInterface.create_basement_bottom(hu_raster_path, surface_raster_file, basement_top_raster_file, basement_bot_raster_file, output_path, raster_driver=raster_driver)

    def read_polyline(self, filename, path=None):

        return self.GISInterface.read_polyline(filename, path)

    def read_points_data(self, filename, path=None):

        return self.GISInterface.read_points_data(filename, path)

    def _getXYpairs(self, points_obj, feature_id=None):

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
        layers = len(rasters) / 2
        points_layer = np.full((layers, len_list), False, dtype=bool)

        return self.GISInterface.map_points_to_raster_layers(points, depths, rasters, points_layer)


    def add2register(self, addition):

        self.model_register += addition

    def writeRegister2file(self):

        with open(self.out_data_folder + 'model_register.dat') as f:
            for item in self.model_register:
                f.write(item)

    def points2shapefile(self, points_array, shapefile_name):
        """ Needs writing ... """
        self.GISInterface.points2shapefile(points_array, shapefile_name)

    def package_data(self):
        """ 
        Option to save all important attributes of DataBuilder class to 
        allow quick loading of data that may have required transforms and 
        processing in its orignial state.

        This will include GIS type objects ... maybe
        """
        target_attr = self.target_attr
        packaged_model = {k: self.__dict__[k] for k in self.__dict__ if k in target_attr}

        # Hack to fix up model boundary which contains a gdal object as well:
        packaged_model['model_boundary'] = packaged_model['model_boundary'][0:4]

        self.save_obj(packaged_model, self.out_data_folder_grid + self.name + '_packaged')

    ###########################################################################

    def __exit__(self):
        self.writeRegister2file()
        self.save_mapped_dictionaries()


class ModelBuilderType(object):

    """
    This class contains all the types that are allowed for different
    class attributes in ModelBuilder.
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
