
class ModelInterface(object):
    """
    The ModelInterface is a class to link the GWModelBuilder with
    different modules for handling building and running of specific
    groundater models, e.g. MODFLOW and HGS
    """    
    
    def __init__(self, model_type=None):
        
        self.model_type = model_type
    
        """        
        :param model_type: Type of model being used.
        """