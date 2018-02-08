import cPickle as pickle
import warnings

from HydroModelBuilder.GWModelBuilder import GWModelBuilder


class GWModelManager(object):

    """Class in which groundwater model lives including functions to build and run the model."""

    def __init__(self, model_directory=None):
        self.model_directory = model_directory
        self.default = 'default'
        self.model_register = None
        self.GW_build = {}

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

    # Save and load utility using pickle
    def save_obj(self, obj, filename):
        """Save object to pickle file.

        :param obj: object, to save.
        :param filename: str, filename to save object to.
        """

        filename = filename + ".pkl" if not filename.endswith(".pkl") else filename
        with open(filename, 'wb') as f:
            pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
        # End with
    # End save_obj()

    def load_obj(self, filename):
        """Load object from pickle file.

        :param filename: str, name of file to load object from.
        """
        if filename.endswith('pkl'):
            with open(filename, 'rb') as f:
                return pickle.load(f)
        else:
            raise TypeError('File type not recognised as "pkl"')
        # End if

    def build_GW_model(self, name=None):
        """Generate a groundwater model using `GWModelBuilder`

        :param name: str, name of model. (Default value = None).
                     If `None` given then name shall be 'default' and its number in the model register.
        """

        if name is None:
            name = 'default{}'.format(self.models + 1)
        self.name = name

        self.GW_build[name] = GWModelBuilder(name=None, model_type=None,
                                             mesh_type=None, units=None,
                                             data_folder=None,
                                             out_data_folder=None,
                                             GISInterface=None,
                                             data_format='binary',
                                             target_attr=self.target_attr)
    # End build_GW_model()

    def emulate_GW_model(self, emulation_method):
        """
        :param emulation_method:
        """
        raise NotImplementedError("Not yet implemented.")
        # Create new GW_model using emulation of existing model
        emulation_methods = ['polynomial chaos expansions', 'PCE', 'multi-fidelity', 'M-F']
        if emulation_method in emulation_methods:
            if emulation_method in ['polynomial chaos expansions', 'PCE']:
                pass
            elif emulation_method in ['multi-fidelity', 'M-F']:
                pass
            # End if
        # End if
    # End if

    def updateParameters(self, model_name, parameters):
        """Update parameters in the groundwater model.

        :param model_name: str, name of model to update.
        :param parameters: str or dict, parameter(s) to update. If string, assume it is a file to load data from.
        """
        warnings.warn("DEPRECATED. Use `update_parameters` instead.", DeprecationWarning)
        self.update_parameters(model_name, parameters)
    # End updateParameters()

    def update_parameters(self, model_name, parameters):
        """Update parameters in the groundwater model.

        :param model_name: str, name of model to update.
        :param parameters: str or dict, parameter(s) to update. If string, assume it is a file to load data from.
        """
        params_new = {}
        if type(parameters) == str:
            # Assume string is a filename to load in the format:
            # Header: PARNAME PARVAL
            # Rows:   p1      val
            with open(parameters, 'r') as f:
                text = f.readlines()
                for line in text:
                    par, val = line.strip('\n').split('\t')
                    params_new[par] = val
                # End for
            # End with
        elif type(parameters) == dict:
            # Assume parameter names match those defined in the model
            # e.g. if params are 'p1' and 'p2' expect e.g. {'p1':0.1, 'p2':3.2}
            params_new = parameters
        # End if

        # Check parameter passed against existing parameters in the model
        original_param_names = self.GW_build[model_name].parameters.param.keys()
        differences = list(set(original_param_names) - set(params_new.keys()))
        if len(differences) > 0:
            print 'The following parameters were not matched: ', differences
        # End if

        for key in params_new:
            if key not in original_param_names:
                print('Parameter name passed not matched, hence not assigned: ', key)
            elif key in original_param_names:
                self.GW_build[model_name].parameters.param[key] = params_new[key]
            # End if
        # End for
    # End update_parameters()

    def setupPEST(self, model_name, directory=None, csv_copy=False, excel_copy=False, models_ID=None):
        """
        :param model_name: str,
        :param directory:  (Default value = None)
        :param csv_copy: (Default value = False)
        :param excel_copy: (Default value = False)
        :param models_ID: (Default value = None)
        """
        from HydroModelBuilder.Utilities.PESTInterface.PESTInterface import PESTInterface
        name = self.GW_build[model_name].name
        if not directory:
            directory = self.GW_build[model_name].out_data_folder_grid
        params = self.GW_build[model_name].parameters.param
        obs = self.GW_build[model_name].observations.obs
        obs_grp = self.GW_build[model_name].observations.obs_group
        print "Model for PEST: ", self.GW_build[model_name]
        self.PEST = PESTInterface(name=name, directory=directory, csv_copy=csv_copy,
                                  excel_copy=excel_copy, params=params, obs=obs, obs_grp=obs_grp, models_ID=models_ID)
    # End setupPEST()

    @property
    def models(self):
        """Get the number of associated models"""
        return len(self.GW_build)
    # End models()

    def load_GW_model(self, GW_model, out_data_folder=None):
        """Load groundwater model from pickle file.

        :param GW_model: str, filename of pickled object to load
        :param out_data_folder: str, where outputs are written to. (Default value = None)
        """

        packaged_model = self.load_obj(GW_model)
        pkg_model_name = packaged_model['name']

        self.GW_build[pkg_model_name] = GWModelBuilder(name=packaged_model['name'],
                                                       model_type=packaged_model['model_type'],
                                                       mesh_type=packaged_model['_mesh_type'],
                                                       units=packaged_model['_units'],
                                                       data_folder=packaged_model['data_folder'],
                                                       out_data_folder=packaged_model['out_data_folder'],
                                                       model_data_folder=packaged_model['model_data_folder'],
                                                       GISInterface=None,
                                                       data_format=packaged_model['data_format'],
                                                       target_attr=self.target_attr)

        ref_pkg_model = self.GW_build[pkg_model_name]
        for key in self.target_attr:
            setattr(ref_pkg_model, key, packaged_model[key])
        # End for

        ref_pkg_model.update_meshgen()
    # End load_GW_model()

    def load_GW_models(self, GW_models):
        """Load given list of groundwater model pickle files.

        :param GW_models: list[str], pickle files to load in.
        """

        for model in GW_models:
            self.load_GW_model(model)
        # End for
    # End load_GW_models()
