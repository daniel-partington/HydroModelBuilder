import cPickle as pickle
import os
import shutil
import sys

import numpy as np
import pandas as pd


class ModelInterface(object):
    """
    The ModelInterface is a class to link the GWModelBuilder with
    different modules for handling building and running of specific
    groundater models, e.g. MODFLOW and HGS
    """

    def __init__(self, model_type, types, data_format):
        """
        :param model_type: str, Type of model being used.
        """
        self.model_type = model_type
        self.types = types
        self.data_format = data_format
        self.units = {}
    # End __init__()

    def set_units(self, length='m', time='s', mass='kg'):
        """
        Sets the units for use inside the model.

        :param length: str, length unit
        :param time: str, time unit
        :param mass: str, mass unit
        """
        if length in self.types.length:
            self.units['length'] = length
        else:
            print 'Length unit not recognised, please use one of: ', self.types.length
            print 'Default unit "m" is set'
        # End if

        if time in self.types.time:
            self.units['time'] = time
        else:
            print 'Time unit not recognised, please use one of: ', self.types.time
            print 'Default unit "s" is set'
        # End if

        if mass in self.types.mass:
            self.units['mass'] = mass
        else:
            print 'Mass unit not recognised, please use one of: ', self.types.mass
            print 'Default unit "kg" is set'
        # End if
    # End set_units()

    def check_for_existing(self, fn):
        """
        Function to determine if input files have previously been processed
        and if so to do nothing unless flagged otherwise. This is done by
        checking the output data path.

        :param fn:
        """
        filename_suffixes = ['_model', '_grid']
        for suffix in filename_suffixes:
            if os.path.isfile(os.path.join(self.out_data_folder, fn[:-4] +
                                           suffix + fn[-4:])):
                print 'found processed file'
            # End if
        # End for
    # End check_for_existing()

    def save_array(self, filename, array):
        if self.data_format == self.types.data_formats[0]:  # if it is 'ascii'
            np.savetxt(filename, array)
        elif self.data_format == self.types.data_formats[1]:  # if it is 'binary'
            np.save(filename, array)
        else:
            print 'Data format not recognised, use "binary" or "ascii"'
        # End if
    # End save_array()

    def load_array(self, array_file):
        if array_file.endswith('txt'):
            return np.loadtxt(array_file)
        elif array_file.endswith('npy') or array_file.endswith('npz'):
            return np.load(array_file)
        else:
            raise TypeError('File type not recognised as "txt", "npy" or "npz" \n {}'.format(array_file))
        # End if
    # End load_array()

    def save_dataframe(self, filename, df):
        df.to_hdf(filename + '.h5', 'table')
    # End save_dataframe()

    def load_dataframe(self, filename):
        if filename.endswith('.h5'):
            return pd.read_hdf(filename, 'table')
        else:
            raise TypeError('File type not recognised as "h5": {}'.format(filename))
        # End if
    # End load_dataframe()

    def save_obj(self, obj, filename):
        with open(filename + '.pkl', 'wb') as f:
            pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
        # End with
    # End save_obj()

    def load_obj(self, filename):
        if filename.endswith('.pkl'):
            with open(filename, 'rb') as f:
                print "Loading: ", f, filename
                p = pickle.load(f)
                return p
            # End with
        else:
            raise TypeError('File type not recognised as "pkl": {}'.format(filename))
        # End if
    # End load_obj()

    def flush(self, mode=None):
        if mode == 'data':
            folder = self.out_data_folder
        elif mode == 'model':
            folder = self.model_data_folder
        else:
            raise ValueError('Expected mode to be either "data" or "model" but got: {}'.format(mode))
        # End if

        if folder == None:
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
        """
        Option to save all important attributes of model class to allow quick
        loading and manipulation of model without separate data and grid generation
        building scripts

        This will not include any GIS type objects
        """
        packaged_model = {k: attrs[k] for k in attrs if k in target_attrs}

        # Hack to fix up model boundary which contains a gdal object as well:
        packaged_model['_model_boundary'] = packaged_model['_model_boundary'][0:4]

        self.save_obj(packaged_model, os.path.join(out_dir, name))
    # End package_model()

# End ModelInterface()
