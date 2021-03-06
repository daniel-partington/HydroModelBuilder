import cPickle as pickle
import os
import shutil

import numpy as np
import pandas as pd


class ModelInterface(object):
    """The ModelInterface is a class to link the GWModelBuilder with
    different modules for handling building and running of specific
    groundater models, e.g. MODFLOW and HGS"""

    def __init__(self, model_type, types, data_format):
        """
        :param model_type: str, type of model being used.
        :param types: object, ModelBuilderType
        :param data_format: str, one of 'ascii' or 'binary'
        """
        assert data_format in ['ascii', 'binary'], "Unknown data format '{}'".format(data_format)
        self.model_type = model_type
        self.types = types
        self.data_format = data_format
        self.units = {}
    # End __init__()

    def set_units(self, length='m', time='s', mass='kg'):
        """Sets the units of measurement used inside the model.

        :param length: str, length unit. Defaults to meters. (Default value = 'm')
        :param time: str, time unit. Defaults to seconds. (Default value = 's')
        :param mass: str, length unit. (Default value = 'kg')
        """
        assert length in self.types.length, "Length unit not recognised, use one of {}".format(self.types.length)
        self.units['length'] = length

        assert time in self.types.time, "Time unit not recognised, use one of {}".format(self.types.time)
        self.units['time'] = time

        assert mass in self.types.mass, "Mass unit not recognised, use one of {}".format(self.types.mass)
        self.units['time'] = time
    # End set_units()

    def check_for_existing(self, fn):
        """Determine if previously processed input files exist in the folder indicated by the `out_data_folder` property.

        :param fn: str, filename to use including extension

        :returns: None, just prints out whether a processed file was found.
        """

        file_ext_split = os.path.splitext(fn)
        filename_suffixes = ['_model', '_grid']
        for suffix in filename_suffixes:
            if os.path.isfile(os.path.join(self.out_data_folder, file_ext_split[0] +
                                           suffix + file_ext_split[1])):
                print('Found processed file')
            # End if
        # End for
    # End check_for_existing()

    def save_array(self, filename, array):
        """Save array data.

        :param filename: str, name of file to save data to
        :param array: numpy array, data to save to file
        """

        if self.data_format in self.types.data_formats:
            # Extract and run associated function that saves data for this format
            save_func = self.types.data_formats[self.data_format]
            save_func(filename, array)
        else:
            print("Data format not recognised, use one of {}".format(self.types.data_formats.keys()))
        # End if
    # End save_array()

    def load_array(self, array_file):
        """Load array data from file.

        :param array_file: str, name of file to load data from, including path and extension
        """

        if array_file.endswith('txt'):
            func = np.loadtxt
        elif array_file.endswith('npy') or array_file.endswith('npz'):
            func = np.load
        else:
            raise TypeError('File type not recognised as "txt", "npy" or "npz" \n {}'.format(array_file))
        # End if

        return func(array_file)
    # End load_array()

    def save_dataframe(self, filename, df):
        """
        :param filename:

        :param df:
        """

        df.to_hdf(filename + '.h5', 'table')
    # End save_dataframe()

    def load_dataframe(self, filename):
        """
        :param filename:
        """

        if filename.endswith('.h5'):
            return pd.read_hdf(filename, 'table')
        else:
            raise TypeError('File type not recognised as "h5": {}'.format(filename))
        # End if
    # End load_dataframe()

    def save_obj(self, obj, filename):
        """
        :param obj:
        :param filename:
        """

        with open(filename + '.pkl', 'wb') as f:
            pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
        # End with
    # End save_obj()

    def load_obj(self, filename):
        """
        :param filename:
        """

        if filename.endswith('.pkl'):
            with open(filename, 'rb') as f:
                print "Loading: ", filename
                p = pickle.load(f)
                return p
            # End with
        else:
            raise TypeError('File type not recognised as "pkl": {}'.format(filename))
        # End if
    # End load_obj()

    def flush(self, mode=None):
        """
        :param mode:  (Default value = None)
        """

        if mode == 'data':
            folder = self.out_data_folder
        elif mode == 'model':
            folder = self.model_data_folder
        else:
            raise ValueError('Expected mode to be either "data" or "model" but got: {}'.format(mode))
        # End if

        if not folder:
            raise ValueError('No folder set, so no flushing')
        # End if
        for the_file in os.listdir(folder):
            file_path = os.path.join(folder, the_file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
                # End if
            except Exception as e:
                print(e)
            # End try
        # End for
    # End flush()

    def package_model(self, out_dir, name, attrs, target_attrs):
        """Option to save all important attributes of model class to allow quick
        loading and manipulation of model without separate data and grid generation
        building scripts

        This will not include any GIS type objects

        :param out_dir:

        :param name:

        :param attrs:

        :param target_attrs:
        """

        packaged_model = {k: attrs[k] for k in attrs if k in target_attrs}

        # Hack to fix up model boundary which contains a gdal object as well:
        packaged_model['_model_boundary'] = packaged_model['_model_boundary'][0:4]
        self.save_obj(packaged_model, os.path.join(out_dir, name))
    # End package_model()

# End ModelInterface()
