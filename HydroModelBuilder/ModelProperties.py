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
    # End __init__()

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
    # End set_temporal_components()
# End ModelTime()


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
        self.bc_types = ['river', 'river_flow', 'wells', 'recharge', 'rainfall',
                         'head', 'drain', 'channel', 'general head']
    # End __init__()

    def create_model_boundary_condition(self, bc_name, bc_type, bc_static=True, bc_parameter=None):
        if bc_type not in self.bc_types:
            raise TypeError('bc_type not recognised, use one of: {}'.format(self.bc_types))
        # End if

        self.bc[bc_name] = {}
        self.bc[bc_name]['bc_type'] = bc_type
        self.bc[bc_name]['bc_static'] = bc_static
    # End create_model_boundary_condition()

    def assign_boundary_array(self, bc_name, bc_array):
        if bc_name in self.bc:
            self.bc[bc_name]['bc_array'] = bc_array
        else:
            sys.exit('No boundary condition with name: {}'.format(bc_name))
        # End if
    # End assign_boundary_array()
# End ModelBoundaries()


class ModelProperties(object):
    """
    This class is used to set all of the parameters which can then be easily
    accessed for modification for model pertubations
    """

    def __init__(self):
        self.properties = {}
        self.prop_types = ['Kh', 'Kv', 'SS', 'Sy']
    # End __init__()

    def assign_model_properties(self, prop_type, value):
        if prop_type in self.prop_types:
            self.properties[prop_type] = value
        else:
            print("{} not in {}".format(prop_type, self.prop_types))
            sys.exit('Property type not recognised')
        # End if
    # End assign_model_properties()
# End ModelProperties()


class ModelParameters(object):
    """
    This class is used to set all of the parameters which can then be easily
    accessed for modification for model pertubations
    """

    def __init__(self):
        self.param = {}
        self.param_set = {}
    # End __init__()

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
    # End create_model_parameter()

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
    # End parameter_options()

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

        for i in xrange(num_parameters):
            name_i = name + str(i)
            self.param[name_i] = {}
            # Test if value is passed as single value or list ... would be nice to accept np arrays here later
            if type(value) == list:
                self.param[name_i]['PARVAL1'] = value[i]
            else:
                self.param[name_i]['PARVAL1'] = value
            # end if
            if i == 0:
                self.param_set[name] = [name_i]
            else:
                self.param_set[name] += [name_i]
            # End if
        # End for
    # End create_model_parameter_set()

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
        # End for
    # End parameter_options_set()
# End ModelParameters()


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
    # End __init__()

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
    # End set_as_observations()

    def collate_observations(self):
        for name in self.obs_group.keys():
            ts = self.obs_group[name]['time_series']
            ts['obs_map'] = 'null'
            for ob in ts.iterrows():
                if ob[1]['active'] == True:
                    self.obs['ob' + str(self.obID)] = ob[1]['value']
                    ts.set_value(ob[0], 'obs_map', 'ob' + str(self.obID))
                    self.obID += 1
                # End if
            # End for
        # End for
    # End collate_observations()

    def check_obs(self):
        obs_nodes = []
        for ob in self.observations.obs.keys():
            obs_nodes += self.observations.obs[ob]
        # End for
    # End check_obs()
# End ModelObservations()


class ModelInitialConditions(object):
    """
    This class is used to store all of the initial conditions for different model
    domains
    """

    def __init__(self):
        self.ic_data = {}
    # End __init__()

    def set_as_initial_condition(self, name, ic_data):
        self.ic_data[name] = ic_data
        return self.ic_data
    # End set_as_initial_condition()
# End ModelInitialConditions


class ModelFeature(object):
    """
    This class defines typical features that might be represented in
    a GW model, the definition of which is assigned to particular arrays
    which can then be passed to the appropriate hydrological model.
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
    # End __init__()
# End ModelFeature()


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
    # End __init__()
# End ModelBuilderType()


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
    # End __init__()

    def SetModflowArrays(self):
        self.layer_order = self.layering_orders[0]
        self.array_order = self.array_ordering[0]
    # End SetModflowArrays()

    def SetHGSArrays(self):
        self.layer_order = self.layering_orders[1]
        self.array_order = self.array_ordering[1]
    # End SetHGSArrays()
# End ArrayOrdering()
