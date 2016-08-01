import sys
import pickle

from HydroModelBuilder.GW_module import GWModelBuilder

class GWModelManager(object):
    
    """ 
    Class in which groundwater model lives including functions to build and
    run the model.
    
    """

    def __init__(self, model_directory=None):
        self.model_directory = model_directory
        self.models = 0
        self.default = 'default'
        self.model_register = None
        self.GW_build = {}

    # Save and load utility using pickle
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
        #end if

    def build_model_database(self):
        pass
    
    def build_GW_model(self, name=None):
        if name == None:
            name ='default' + self.models
        self.name = name        
        self.models += 1 
        self.GW_build[self.models] = GWModelBuilder(name=None, model_type=None, 
                                                    mesh_type=None, units=None, 
                                                    data_folder=None, 
                                                    out_data_folder=None, 
                                                    GISInterface=None, 
                                                    data_format='binary')

    def emulate_GW_model(self, emulation_method):
        # Create new GW_model using emulation of existing model    
        emulation_methods = ['polynomial chaos expansions', 'PCE', 'multi-fidelity', 'M-F']
        if emulation_method in emulation_methods:
            if emulation_method in ['polynomical chaos expansions', 'PCE']:
                pass
            elif emulation_method in ['multi-fidelity', 'M-F']:
                pass

    def load_GW_model(self, GW_model, out_data_folder=None):
        packaged_model = self.load_obj(GW_model)

        self.models += 1 
        self.GW_build[self.models] = GWModelBuilder(name=packaged_model['name'], 
                                                    model_type=packaged_model['model_type'], 
                                                    mesh_type=packaged_model['mesh_type'], 
                                                    units=packaged_model['units'],
                                                    data_folder=packaged_model['data_folder'], 
                                                    out_data_folder=packaged_model['out_data_folder'], 
                                                    GISInterface=None, 
                                                    data_format=packaged_model['data_format'])

        # Rename model in dictionary by it's model builder name
        self.GW_build[packaged_model['name']] = self.GW_build.pop(self.models)
        # Load in all of the objects except for GIS objects to the model builder class        
        
        #self.GW_build[self.models].name = packaged_model['name']
        #self.GW_build[self.models].model_type = packaged_model['model_type']
        #self.GW_build[self.models].mesh_type = packaged_model['mesh_type']
        #self.GW_build[self.models].data_format = packaged_model['data_format'] 
        #self.GW_build[self.models].out_data_folder = packaged_model['out_data_folder']
        #self.GW_build[self.models].data_folder = packaged_model['data_folder']
        #self.GW_build[self.models].out_data_folder = packaged_model['out_data_folder']
        
        self.GW_build[packaged_model['name']].array_ordering = packaged_model['array_ordering']
        self.GW_build[packaged_model['name']].boundaries = packaged_model['boundaries']
        self.GW_build[packaged_model['name']].parameters = packaged_model['parameters']
        self.GW_build[packaged_model['name']].observations = packaged_model['observations']
        self.GW_build[packaged_model['name']].initial_conditions = packaged_model['initial_conditions']
        self.GW_build[packaged_model['name']].out_data_folder_grid = packaged_model['out_data_folder_grid']
        self.GW_build[packaged_model['name']].model_boundary = packaged_model['model_boundary']
        self.GW_build[packaged_model['name']].boundary_poly_file = packaged_model['boundary_poly_file']
        #self.GW_build[self.models].data_boundary = packaged_model['data_boundary']
        self.GW_build[packaged_model['name']].boundary_data_file = packaged_model['boundary_data_file']
        #self.GW_build[self.models].model_mesh = packaged_model['model_mesh']
        self.GW_build[packaged_model['name']].model_mesh_centroids = packaged_model['model_mesh_centroids']
        self.GW_build[packaged_model['name']].model_mesh3D = packaged_model['model_mesh3D']
        self.GW_build[packaged_model['name']].model_mesh3D_centroids = packaged_model['model_mesh3D_centroids']
        self.GW_build[packaged_model['name']].model_layers = packaged_model['model_layers']
        self.GW_build[packaged_model['name']].model_features = packaged_model['model_features']
        self.GW_build[packaged_model['name']].polyline_mapped = packaged_model['polyline_mapped']
        self.GW_build[packaged_model['name']].points_mapped = packaged_model['points_mapped']
        self.GW_build[packaged_model['name']].model_register = packaged_model['model_register']
        self.GW_build[packaged_model['name']].base_data_register = packaged_model['base_data_register']
        self.GW_build[packaged_model['name']].gridded_data_register = packaged_model['gridded_data_register']
                              
        # Temporary package variables that need embedding elsewhere ...
        self.GW_build[packaged_model['name']].gridHeight = packaged_model['gridHeight']
        self.GW_build[packaged_model['name']].gridWidth = packaged_model['gridWidth']
        
        
    def load_GW_models(self, GW_models):
        for model in GW_models:
            self.model = self.load_GW_model(model)
    
    def GW_model_runner(self):
        pass
    
    def create_GW_model_copy(self):
        pass