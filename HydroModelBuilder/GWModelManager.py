import cPickle as pickle

from HydroModelBuilder.GWModelBuilder import GWModelBuilder


class GWModelManager(object):

    """
    Class in which groundwater model lives including functions to build and
    run the model.
    """

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
        with open(filename + '.pkl', 'wb') as f:
            pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
        # End with
    # End save_obj()

    def load_obj(self, filename):
        if filename.endswith('pkl'):
            with open(filename, 'rb') as f:
                return pickle.load(f)
        else:
            raise TypeError('File type not recognised as "pkl"')
        # End if

    def build_model_database(self):
        pass

    def build_GW_model(self, name=None):
        if name is None:
            name = 'default{}'.format(self.models)
        self.name = name

        self.GW_build[name] = GWModelBuilder(name=None, model_type=None,
                                             mesh_type=None, units=None,
                                             data_folder=None,
                                             out_data_folder=None,
                                             GISInterface=None,
                                             data_format='binary',
                                             target_attr=self.target_attr)

    def emulate_GW_model(self, emulation_method):
        # Create new GW_model using emulation of existing model
        emulation_methods = ['polynomial chaos expansions', 'PCE', 'multi-fidelity', 'M-F']
        if emulation_method in emulation_methods:
            if emulation_method in ['polynomial chaos expansions', 'PCE']:
                pass
            elif emulation_method in ['multi-fidelity', 'M-F']:
                pass

    def updateParameters(self, model_name, parameters):

        params_new = {}
        if type(parameters) == str:
            """
            Assume string is a filename to load in the format:
            Header: PARNAME PARVAL
            Rows:   p1      val
            """
            with open(parameters, 'r') as f:
                text = f.readlines()
                for line in text:
                    par, val = line.strip('\n').split('\t')
                    params_new[par] = val

        elif type(parameters) == dict:
            """
            Assume parameter names match those defined in the model
            e.g. if params are 'p1' and 'p2' expect e.g. {'p1':0.1, 'p2':3.2}
            """
            params_new = parameters

        # Check parameter passed against existing parameters in the model
        original_param_names = self.GW_build[model_name].parameters.param.keys()
        differences = list(set(original_param_names) - set(params_new.keys()))
        if len(differences) > 0:
            print 'The following parameters were no matched: ', differences
        # End if
        for key in params_new.keys():
            if key not in original_param_names:
                print 'Parameter name passed not matched, hence not assigned: ', key
            elif key in original_param_names:
                self.GW_build[model_name].parameters.param[key] = params_new[key]
            # End if
        # end for

    def setupPEST(self, model_name, directory=None, csv_copy=False, excel_copy=False, models_ID=None):
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
        """Load groundwater model."""
        packaged_model = self.load_obj(GW_model)

        self.GW_build[packaged_model['name']] = GWModelBuilder(name=packaged_model['name'],
                                                               model_type=packaged_model['model_type'],
                                                               mesh_type=packaged_model['_mesh_type'],
                                                               units=packaged_model['_units'],
                                                               data_folder=packaged_model['data_folder'],
                                                               out_data_folder=packaged_model['out_data_folder'],
                                                               model_data_folder=packaged_model['model_data_folder'],
                                                               GISInterface=None,
                                                               data_format=packaged_model['data_format'],
                                                               target_attr=self.target_attr)

        ref_pkg_model = self.GW_build[packaged_model['name']]
        for key in self.target_attr:
            setattr(ref_pkg_model, key, packaged_model[key])
        # End for

        ref_pkg_model.update_meshgen()
    # End load_GW_model()

    def load_GW_models(self, GW_models):
        for model in GW_models:
            self.model = self.load_GW_model(model)

    def GW_model_runner(self):
        pass

    def create_GW_model_copy(self):
        pass
