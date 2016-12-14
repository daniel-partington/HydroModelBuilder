import os

import numpy as np
import matplotlib.pyplot as plt

import flopy
import flopy.utils.binaryfile as bf

class ModflowModel(object):

    def __init__(self, model_data, data_folder=None, **kwargs):

        """
        :param name: model_data: Object containing all the data for the model
        """

        self.model_data = model_data
        self.name = model_data.name
        if data_folder == None:
            self.data_folder = os.getcwd() + os.path.sep + 'model_' + self.name + os.path.sep
        else:
            self.data_folder = data_folder + os.path.sep + 'model_' + self.name + os.path.sep
        #end if
            
        self.executable = r".\MODFLOW-NWT_64.exe"

        self.nlay = self.model_data.model_mesh3D[0].shape[0]-1
        self.nrow = self.model_data.model_mesh3D[0].shape[1]
        self.ncol = self.model_data.model_mesh3D[0].shape[2]
        self.delr = self.model_data.gridHeight
        self.delc = self.model_data.gridWidth
        self.top = self.model_data.model_mesh3D[0][0]
        self.botm = self.model_data.model_mesh3D[0][1:] 
        self.xul = self.model_data.model_boundary[0]
        self.yul = self.model_data.model_boundary[3]
        
        if self.model_data.model_time.t['steady_state'] == True:
            self.nper = 1
            self.perlen = 40000 * 365 #8260000#65 #000000 
            self.nstp = 1 #0
            self.steady = True #False # False #True #False#True #False
            #self.start_datetime = self.model_data.model_time.t['start_time']
        else:
            self.nper = self.model_data.model_time.t['steps']
            self.perlen = [x.total_seconds()/86400. for x in self.model_data.model_time.t['intervals']] #365000 
            self.nstp = 1 #10
            self.steady = False 
            self.start_datetime = self.model_data.model_time.t['start_time']
        
        # Initial data:
        self.strt = self.model_data.initial_conditions.ic_data["Head"] #self.model_data.model_mesh3D[0][1:] + 20.         
        
        self.hk = self.model_data.properties.properties['Kh']
        self.hk[self.hk == -1] = 1.
        self.vka = self.model_data.properties.properties['Kv']
        self.sy = self.model_data.properties.properties['Sy']
        self.ss = self.model_data.properties.properties['SS']

        #Set all other kwargs as class attributes
        for key, value in kwargs.items():
            setattr(self, key, value)
        #End For

        #self.flowpy_params['model_dir'] = os.path.join(data_folder, "MF_IO", name)

    #End init()

    def createDiscretisation(self):

        # Create discretisation object
        self.dis = flopy.modflow.ModflowDis(self.mf, 
                                            nlay = self.nlay,
                                            nrow = self.nrow,
                                            ncol = self.ncol,
                                            delr = self.delr,
                                            delc = self.delc,
                                            top = self.top,
                                            botm = self.botm,
                                            xul = self.xul,
                                            yul = self.yul,
                                            nper = self.nper,
                                            perlen = self.perlen,
                                            nstp = self.nstp,
                                            steady = self.steady) #,
                                            #start_datetime = self.start_datetime)
    #End createDiscretisation()

    def setupBASPackage(self, nlay, nrow, ncol):
    
        # Variables for the BAS package
        #ibound = np.nan*np.empty((nlay, nrow, ncol), dtype=np.int32)
        ibound = self.model_data.model_mesh3D[1]        
        ibound[ibound == -1] = 0
        #ibound[ibound > 0] = 1
        
        strt = self.strt #np.nan*np.empty((nlay, nrow, ncol), dtype=np.float32)

        self.bas = flopy.modflow.ModflowBas(self.mf, ibound=ibound, strt=strt)

    #End setupBASPackage()

    def setupNWTpackage(self):
        # Specify NWT settings
        self.nwt = flopy.modflow.ModflowNwt(self.mf, headtol=1E-4, fluxtol=1.0E1, linmeth=2, thickfact=1E-7, maxiterout=100000, options='COMPLEX')

    #end setupNWTpackage

    def setupPCGpackage(self):
        self.pcg = flopy.modflow.ModflowPcg(self.mf)  # if using mf 2005
        
    # end setupPCGpachage    

    def setupUPWpackage(self ):
        """
        This function is used to set up the upstream weighting package for the
        NWT solver in MODFLOW
        """        
        # Add UPW package to the MODFLOW model to represent aquifers
        self.upw = flopy.modflow.ModflowUpw(self.mf, 
                                            hk=self.hk, 
                                            vka=self.vka, 
                                            sy=self.sy, 
                                            ss=self.ss, 
                                            hdry=-999.9, 
                                            laytyp=1)

    # end setupUPWpackage
    
    def createRIVpackage(self, lrcd=None):
        self.riv = flopy.modflow.ModflowRiv(self.mf, ipakcb=53, stress_period_data=lrcd)

    # end createRIVpackage        
       
    def createSTRpackage(self, STRvariables=None):
        self.str = flopy.modflow.ModflowStr(self.mf, ipakcb=53, stress_period_data=STRvariables)

    def createDRNpackage(self, lrcsc=None):
        self.drn = flopy.modflow.ModflowDrn(self.mf, ipakcb=53, stress_period_data=lrcsc)    
    
    def createGHBpackage(self, lrcsc=None):
        self.ghb = flopy.modflow.ModflowGhb(self.mf, ipakcb=53, stress_period_data=lrcsc)
       
    def createRCHpackage(self, rchrate=None):
        # Add RCH package to the MODFLOW model to represent recharge
        # rchrate  = 1.0E-3 * np.random.rand(self.nrow, self.ncol)
        self.rch = flopy.modflow.ModflowRch(self.mf, ipakcb=53, rech=rchrate, nrchop=3)
    #end createRCHpackage
        
    def createWELpackage(self, lrcq=None):
        # Add WEL package to the MODFLOW model to represent pumping wells
        # Expects a dictionary with an array of well location lay row col and flux
        # at each stress period, e.g.        
        #lrcq = {}
        #lrcq[0] = [[0, 7, 7, -100.]] # layer, row, column, flux
        self.wel = flopy.modflow.ModflowWel(self.mf, ipakcb=53, stress_period_data=lrcq)

    # end createWElpackage    
    
    def createOCpackage(self):
        # Add OC package to the MODFLOW model
        spd = {(0, 0): ['save head', 'print budget', 'save budget'], (0, 1): ['save head', 'print budget', 'save budget']}
        self.oc = flopy.modflow.ModflowOc(self.mf, stress_period_data=spd, cboufm='(20i5)')

    # end createOCpackage

    def createLMTpackage(self):
        # Add LMT package to the MODFLOW model to allow linking with MT3DMS
        self.lmt = flopy.modflow.ModflowLmt(self.mf)

    # end createLMTpackage
        
    def finaliseModel(self):
        # Write the MODFLOW model input files
        self.mf.write_input()        
    # end finaliseModel

    def buildMODFLOW(self):

        self.mf = flopy.modflow.Modflow(self.name, exe_name=self.executable, 
                                        model_ws=self.data_folder, version='mfnwt')

        self.createDiscretisation()
        self.setupBASPackage(self.nlay, self.nrow, self.ncol)
        self.setupNWTpackage()        
        self.setupUPWpackage()    
        self.createOCpackage()
        
        for boundary in self.model_data.boundaries.bc:
            if self.model_data.boundaries.bc[boundary]['bc_type'] == 'recharge':
                self.createRCHpackage(rchrate=self.model_data.boundaries.bc[boundary]['bc_array'])
                
            #if self.model_data.boundaries.bc[boundary]['bc_type'] == 'river':
            #    self.createRIVpackage(self.model_data.boundaries.bc[boundary]['bc_array'])
                
            #if self.model_data.boundaries.bc[boundary]['bc_type'] == 'wells':
            #    self.createWELpackage(self.model_data.boundaries.bc[boundary]['bc_array'])

            if self.model_data.boundaries.bc[boundary]['bc_type'] == 'drain':
                self.createDRNpackage(self.model_data.boundaries.bc[boundary]['bc_array'])

            if self.model_data.boundaries.bc[boundary]['bc_type'] == 'general head':
                self.createGHBpackage(self.model_data.boundaries.bc[boundary]['bc_array'])

        river_exists = False
        wells_exist = False
        for boundary in self.model_data.boundaries.bc:
            if self.model_data.boundaries.bc[boundary]['bc_type'] == 'river':
                river_exists = True
            elif self.model_data.boundaries.bc[boundary]['bc_type'] == 'wells':
                wells_exist = True
        
        if river_exists:
            river = {}
            river[0] = []
            for boundary in self.model_data.boundaries.bc:
                if self.model_data.boundaries.bc[boundary]['bc_type'] == 'river':
                    time_key = self.model_data.boundaries.bc[boundary]['bc_array'].keys()[0]
                    river[0] += self.model_data.boundaries.bc[boundary]['bc_array'][time_key]
                if self.model_data.boundaries.bc[boundary]['bc_type'] == 'channel':
                    time_key = self.model_data.boundaries.bc[boundary]['bc_array'].keys()[0]
                    river[0] += self.model_data.boundaries.bc[boundary]['bc_array'][time_key]
 
            self.createRIVpackage(river)

        if wells_exist:
            wel = {}
            wel[0] = []
            for boundary in self.model_data.boundaries.bc:
                if self.model_data.boundaries.bc[boundary]['bc_type'] == 'wells':
                    time_key = self.model_data.boundaries.bc[boundary]['bc_array'].keys()[0]
                    wel[0] += self.model_data.boundaries.bc[boundary]['bc_array'][time_key]

            self.createWELpackage(wel)
            
            
        # IF transport then
        self.createLMTpackage()
            
        self.finaliseModel()
    #End buildMODFLOW        

    def checkMODFLOW(self):
        self.mf.check()
    #End checkMODFLOW        

    def runMODFLOW(self):

        success, buff = self.mf.run_model()
        if not success:
            raise Exception('MODFLOW did not terminate normally.')
        #End if
    #End runMODFLOW()
    
    
    #**************************************************************************
    #**************************************************************************
    #**************************************************************************
    #**** CHECKING SUCCESS OF MODEL RUN ***************************************
    #**************************************************************************
    #**************************************************************************
    #**************************************************************************


    def checkCovergence(self, path=None, name=None):
        converge_fail_options = ["****FAILED TO MEET SOLVER CONVERGENCE CRITERIA IN TIME STEP", # Clear statement of model fail in list file 
                                 " PERCENT DISCREPANCY =         200.00", # Convergence but extreme discrepancy in results
                                 " NaN " # Something big went wrong but somehow convergence was reached?
                                 ]
        if path:
            with open(os.path.join(path,name) + '.list', 'r') as f:
                list_file = f.read()
            for converge_fail in converge_fail_options:
                if converge_fail in list_file:
                    print "*** Convergence failure ***"            
                    print os.getcwd()            
                    import datetime                
                    now = datetime.datetime.now().strftime("%I%M%p%B%d%Y")
                    with open(os.path.join(self.data_folder, "converge_fail_%s.txt" %now), 'w') as fail:
                        fail.write("Model did not converge, @ %s" %datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y"))
                    return False
                #end if
            #end for
            return True

        else:
            with open(self.data_folder + self.name + '.list', 'r') as f:
                list_file = f.read()
            for converge_fail in converge_fail_options:
                if converge_fail in list_file:
                    print "*** Convergence failure ***"            
                    print os.getcwd()            
                    import datetime                
                    now = datetime.datetime.now().strftime("%I%M%p%B%d%Y")
                    with open(os.path.join(self.data_folder, "converge_fail_%s.txt" %now), 'w') as fail:
                        fail.write("Model did not converge, @ %s" %datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y"))
                    return False
                #end if
            #end for
            return True


    #**************************************************************************
    #**************************************************************************
    #**************************************************************************
    #**** PROCESSING OF MODEL RUN *********************************************
    #**************************************************************************
    #**************************************************************************
    #**************************************************************************

    def getRiverFlux(self, name):
        cbbobj = self.importCbb()
        riv_flux = cbbobj.get_data(text='RIVER LEAKAGE', full3D=True)            
        times = cbbobj.get_times()
        if type(times) == str:
            str_pers = 1
        else:
            str_pers = len(times)
        if self.model_data.boundaries.bc[name]['bc_type'] == 'river':
            time_key = self.model_data.boundaries.bc[name]['bc_array'].keys()[0]
            river = self.model_data.boundaries.bc[name]['bc_array'][time_key]
        else:
            print 'Not a river boundary'

        riv_exchange = {}
        for per in range(str_pers):
            riv_exchange[per] = []
            for cell in river:
                (l, r, c) = (cell[0], cell[1], cell[2])            
                if str_pers > 1:                
                    riv_exchange[per] += [[riv_flux[per][l][r][c], (l,r,c)]]
                else:
                    riv_exchange[per] += [[riv_flux[0][l][r][c], (l,r,c)]]

        return riv_exchange

    def getRiverFluxNodes(self, nodes):
        cbbobj = self.importCbb()
        riv_flux = cbbobj.get_data(text='RIVER LEAKAGE', full3D=True)            
        times = cbbobj.get_times()
        if type(times) == str:
            str_pers = 1
        else:
            str_pers = len(times)
        #end if
        river = nodes

        riv_exchange = {}
        for per in range(str_pers):
            riv_exchange[per] = []
            for cell in river:
                (l, r, c) = (0, cell[0], cell[1])            
                if str_pers > 1:                
                    riv_exchange[per] += [[riv_flux[per][l][r][c], (l,r,c)]]
                else:
                    riv_exchange[per] += [[riv_flux[0][l][r][c], (l,r,c)]]

                                          
        riv_exchange = np.array([x[0] for x in riv_exchange[0] if type(x[0]) == np.float32]).sum()
                                          
        return riv_exchange

    def getAverageDepthToGW(self, mask=None):
        
        headobj = self.importHeads()
        times = headobj.get_times()
        head = headobj.get_data(totim=times[-1])

        if mask is None:
            return np.mean(self.top - head[0])
        else:
            return np.mean(self.top[mask] - head[0][mask])
        
    def ConcsByZone(self, concs):
        
        self.model_data.model_mesh3D[1]

        concs_zoned = [np.full(self.model_data.model_mesh3D[1].shape[1:3], np.nan)] * int(np.max(self.model_data.model_mesh3D[1]))
        
        for zone in range(int(np.max(self.model_data.model_mesh3D[1]))):
            temp_concs = np.array([np.full(self.model_data.model_mesh3D[1].shape[1:3], np.nan)] * self.model_data.model_mesh3D[1].shape[0])
            # for each layer in mesh get the heads from zone and average:
            for layer in range(self.model_data.model_mesh3D[1].shape[0]):
                temp_concs[layer][self.model_data.model_mesh3D[1][layer] == float(zone+1)] = concs[layer][self.model_data.model_mesh3D[1][layer] == float(zone+1)]
            masked_temp_concs = np.ma.masked_array(temp_concs, np.isnan(temp_concs))
            concs_zoned[zone] = np.mean(masked_temp_concs, axis=0) 
        #end for
        return concs_zoned

    def HeadsByZone(self, heads):
        
        self.model_data.model_mesh3D[1]

        heads_zoned = [np.full(self.model_data.model_mesh3D[1].shape[1:3], np.nan)] * int(np.max(self.model_data.model_mesh3D[1]))
        
        for zone in range(int(np.max(self.model_data.model_mesh3D[1]))):
            temp_heads = np.array([np.full(self.model_data.model_mesh3D[1].shape[1:3], np.nan)] * self.model_data.model_mesh3D[1].shape[0])
            # for each layer in mesh get the heads from zone and average:
            for layer in range(self.model_data.model_mesh3D[1].shape[0]):
                temp_heads[layer][self.model_data.model_mesh3D[1][layer] == float(zone+1)] = heads[layer][self.model_data.model_mesh3D[1][layer] == float(zone+1)]
            masked_temp_heads = np.ma.masked_array(temp_heads, np.isnan(temp_heads))
            heads_zoned[zone] = np.mean(masked_temp_heads, axis=0) 
        #end for
        return heads_zoned

#    def ObservedHeadsByZone(self, obs_set, head):
#        self.obs_heads_by_zone = [] * int(np.max(self.model_data.model_mesh3D[1]))
#        obs_df = self.model_data.observations.obs_group[obs_set]['time_series']
#        obs_df = obs_df[obs_df['active'] == True]
#        sim_map_dict = self.model_data.observations.obs_group[obs_set]['mapped_observations']        
#
#        for observation in obs_df['name']:
#            idx = obs_df[obs_df['name'] == observation].index.tolist()[0]            
#            obs = obs_df.get_value(idx, 'value')
#            self.obs_heads_by_zone[] +=             

    def writeObservations(self):

        # Set model output arrays to None to initialise        
        head = None
        
        with open(self.data_folder + os.path.sep + 'observations_'+ self.model_data.name +'.txt', 'w') as f:
            # Write observation to file         
            for obs_set in self.model_data.observations.obs_group.keys():
                obs_type = self.model_data.observations.obs_group[obs_set]['obs_type']
                # Import the required model outputs for processing            
                if obs_type == 'head':
                    # Check if model outputs have already been imported and if not import                
                    if not head:
                        headobj = self.importHeads()
                        #times = headobj.get_times()
                        head = headobj.get_alldata() #(totim=times[0])
                elif obs_type == 'stage':
                    #stage = self.importStage()
                    pass
                elif obs_type == 'discharge':
                    pass
                #elif obs_type == 'concentration':
                #    # Check if model outputs have already been imported and if not import                
                #   if not head:
                #        concobj = self.importConcs()
                #        conc = concobj.get_alldata() #(totim=times[0])
                #end if
                
                obs_df = self.model_data.observations.obs_group[obs_set]['time_series']
                obs_df = obs_df[obs_df['active'] == True]
                sim_map_dict = self.model_data.observations.obs_group[obs_set]['mapped_observations']        

                for observation in obs_df.index:
                    interval = obs_df['interval'].loc[observation]
                    name = obs_df['name'].loc[observation]
                    x = self.model_data.observations.obs_group[obs_set]['locations']['Easting'].loc[name]            
                    y = self.model_data.observations.obs_group[obs_set]['locations']['Northing'].loc[name]            
                    (x_cell, y_cell) = self.model_data.mesh2centroid2Dindex[(sim_map_dict[name][1], sim_map_dict[name][2])]
                    west = True
                    north = True                    
                    (lay, row, col) = [sim_map_dict[name][0],sim_map_dict[name][1],sim_map_dict[name][2]]                     
                    sim_heads = [head[interval][lay][row][col]]

#                    if x > x_cell:
#                        west = False
#                        west_diff = x-x_cell
#                    if y < y_cell:
#                        north = False
#                        north_diff = y - y_cell
#
#                    sim_zone = self.model_data.model_mesh3D[1][lay][row][col]                     
#                    
#                    if north:     
#                        try:
#                            if self.model_data.model_mesh3D[1][lay][row-1][col] == sim_zone:
#                                sim_heads += [head[interval][lay][row-1][col]]
#                        except IndexError as e:
#                            print e.mmsg
#                        if west:     
#                            try:
#                                if self.model_data.model_mesh3D[1][lay][row-1][col-1] == sim_zone:
#                                    sim_heads += [head[interval][lay][row-1][col-1]]
#                             except IndexError as e:
#                                 print e.mmsg
#                        else:
#                            try:
#                                if self.model_data.model_mesh3D[1][lay][row-1][col+1] == sim_zone:
#                                    sim_heads += [head[interval][lay][row-1][col+1]]
#                            except:
#                                pass
#
#                    else:
#                        try:
#                            if self.model_data.model_mesh3D[1][lay][row+1][col] == sim_zone:
#                                sim_heads += [head[interval][lay][row+1][col]]
#                        except:
#                            pass
#                        if west:     
#                            try:
#                                if self.model_data.model_mesh3D[1][lay][row+1][col-1] == sim_zone:
#                                    sim_heads += [head[interval][lay][row+1][col-1]]
#                            except:
#                                pass
#                        else:
#                            try:
#                                if self.model_data.model_mesh3D[1][lay][row+1][col+1] == sim_zone:
#                                    sim_heads += [head[interval][lay][row+1][col+1]]
#                            except:
#                                pass
#
#
#                    if west:     
#                        try:
#                            if self.model_data.model_mesh3D[1][lay][row][col-1] == sim_zone:
#                                sim_heads += [head[interval][lay][row][col-1]]
#                        except:
#                            pass
#                    else:
#                        try:
#                            if self.model_data.model_mesh3D[1][lay][row][col+1] == sim_zone:
#                                sim_heads += [head[interval][lay][row][col+1]]
#                        except:
#                            pass

                    sim_head = np.mean(sim_heads)
                    f.write('%f\n' %sim_head)

                    
                    
    def getObservation(self, obs, interval, obs_set):
        # Set model output arrays to None to initialise        
        head = None
        
        obs_type = self.model_data.observations.obs_group[obs_set]['obs_type']
        # Import the required model outputs for processing            
        if obs_type == 'head':
            # Check if model outputs have already been imported and if not import                
            if not head:
                headobj = self.importHeads()
                #times = headobj.get_times()
                head = headobj.get_alldata() #(totim=times[0])
        elif obs_type == 'stage':
            #stage = self.importStage()
            pass
        #end if
        sim_map_dict = self.model_data.observations.obs_group[obs_set]['mapped_observations']        
        sim_head = head[interval][sim_map_dict[obs][0]][sim_map_dict[obs][1]][sim_map_dict[obs][2]]
        dtw = self.model_data.model_mesh3D[0][0][sim_map_dict[obs][1]][sim_map_dict[obs][2]] - sim_head        
        return sim_head, dtw
    

    def CompareObservedHead(self, obs_set, simulated, nper=0):
        
        self.comparison = []
        self.comp_zone = []
        self.obs_loc_val_zone = []
        obs_df = self.model_data.observations.obs_group[obs_set]['time_series']
        obs_df = obs_df[obs_df['active'] == True]
        obs_df = obs_df[obs_df['interval'] == nper]
        sim_map_dict = self.model_data.observations.obs_group[obs_set]['mapped_observations']        
        for observation in obs_df['name']: #self.model_data.observations.obs_group[obs_set]['time_series']['name']:
            idx = obs_df[obs_df['name'] == observation].index.tolist()[0]            
            obs = obs_df.get_value(idx, 'value')
            sim = simulated[sim_map_dict[observation][0]][sim_map_dict[observation][1]][sim_map_dict[observation][2]]

            self.comparison += [[obs, sim]] 
            self.comp_zone += self.model_data.model_mesh3D[1][sim_map_dict[observation][0]][sim_map_dict[observation][1]][sim_map_dict[observation][2]]/np.max(self.model_data.model_mesh3D[1])
            
            x = self.model_data.observations.obs_group[obs_set]['locations']['Easting'].loc[observation]            
            y = self.model_data.observations.obs_group[obs_set]['locations']['Northing'].loc[observation]            
            zone = self.model_data.model_mesh3D[1][sim_map_dict[observation][0]][sim_map_dict[observation][1]][sim_map_dict[observation][2]]            
            self.obs_loc_val_zone += [[obs, (x, y), zone, sim]]

    def CompareObserved(self, obs_set, simulated, nper=0):
        
        self.obs_sim_zone = []
        obs_df = self.model_data.observations.obs_group[obs_set]['time_series']
        obs_df = obs_df[obs_df['active'] == True]
        obs_df = obs_df[obs_df['interval'] == nper]
        sim_map_dict = self.model_data.observations.obs_group[obs_set]['mapped_observations']        
        for observation in obs_df['name']: #self.model_data.observations.obs_group[obs_set]['time_series']['name']:
            idx = obs_df[obs_df['name'] == observation].index.tolist()[0]            
            obs = obs_df.get_value(idx, 'value')
            sim = simulated[sim_map_dict[observation][0]][sim_map_dict[observation][1]][sim_map_dict[observation][2]]
            zone = self.model_data.model_mesh3D[1][sim_map_dict[observation][0]][sim_map_dict[observation][1]][sim_map_dict[observation][2]]            
            x = self.model_data.observations.obs_group[obs_set]['locations']['Easting'].loc[observation]            
            y = self.model_data.observations.obs_group[obs_set]['locations']['Northing'].loc[observation]            
            if np.isnan(sim):
                print sim, obs, zone                
                continue
            self.obs_sim_zone += [[obs, sim, zone, x, y]]
        
    
    def importHeads(self, path=None, name=None):
        if path:
            headobj = bf.HeadFile(path + name +'.hds') #, precision='double')
            return headobj
        else:
            self.headobj = bf.HeadFile(self.data_folder + self.name+'.hds')
            return self.headobj

    def importCbb(self):
        self.cbbobj = bf.CellBudgetFile(self.data_folder + self.name+'.cbc')
        return self.cbbobj

    def waterBalance(self):
        import flopy.utils.binaryfile as bf
        import pandas as pd
        
        cbbobj = bf.CellBudgetFile(self.data_folder + self.name+'.cbc')

        water_balance_components = cbbobj.textlist
        water_balance = {}        
        water_balance_summary = {}
        water_bal_summed_titles = []
        water_bal_summed_values = []
        for component in water_balance_components:
            component_stripped = component.lstrip().rstrip()
            if 'FLOW' in component_stripped: continue            
            water_balance[component_stripped] = cbbobj.get_data(text = component, full3D=True)#[-1]            
            #print water_balance[component_stripped]            
            if np.any(water_balance[component_stripped][0] > 0):
                water_balance_summary[component_stripped+'_pos'] = np.sum(water_balance[component_stripped][0][water_balance[component_stripped][0] > 0])
            else:
                water_balance_summary[component_stripped+'_pos'] = 0.0
            if np.any(water_balance[component_stripped][0] < 0):
                water_balance_summary[component_stripped+'_neg'] = np.sum(water_balance[component_stripped][0][water_balance[component_stripped][0] < 0])
            else:
                water_balance_summary[component_stripped+'_neg'] = 0.0
            
            water_bal_summed_titles += component_stripped+'_pos', component_stripped+'_neg'
            water_bal_summed_values += water_balance_summary[component_stripped+'_pos'], water_balance_summary[component_stripped+'_neg']            

        water_bal_summed_titles += ['Error']
        Error = sum(water_bal_summed_values)        
        water_bal_summed_values += [Error]       
        wat_bal_df = pd.DataFrame(water_bal_summed_values, water_bal_summed_titles) # pd.DataFrame.from_dict(water_balance_summary, orient='index')
        wat_bal_df.columns = ['Flux m^3/d']        
        wat_bal_df = wat_bal_df[wat_bal_df['Flux m^3/d'] != 0.]
        print wat_bal_df       

        # Setup params to get water balance aspect ratio looking nice
        #aspect = float(12.5715/((wat_bal_df.max()[0]-wat_bal_df.min()[0])/float(wat_bal_df.shape[1])))        

        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(1, 1, 1) #, aspect=aspect)
        ax.set_title('Water Balance')
        wat_bal_df.plot(kind='bar', ax=plt.gca())        
        ax.grid(True)        
        gridlines = ax.get_xgridlines() #+ ax.get_ygridlines()
        for line in gridlines:
            line.set_linestyle('-')

        fig.subplots_adjust(left=0.1, right=0.9, bottom=0.35, top=0.95, wspace=0.1, hspace=0.12) 
        
    def compareAllObs(self):

        headobj = self.importHeads()
        times = headobj.get_times()        

        scatterx = []        
        scattery = []        
        obs_sim_zone_all = []

        # The definition of obs_sim_zone looks like:
        # self.obs_sim_zone += [[obs, sim, zone, x, y]]
         
        for i in range(self.model_data.model_time.t['steps']):
            head = headobj.get_data(totim=times[i])
            self.CompareObserved('head', head, nper=i)                
            obs_sim_zone_all += self.obs_sim_zone

#        obs_sim_zone_cull = []         
#        for x in obs_sim_zone_all:
#            if abs(x[0]-x[1] > 50):
#                continue
#            else:
#                obs_sim_zone_cull += [x]
#        
#        obs_sim_zone_all = obs_sim_zone_cull

        scatterx = np.array([h[0] for h in obs_sim_zone_all])        
        scattery = np.array([h[1] for h in obs_sim_zone_all])        
        #print np.min(scatterx), np.max(scatterx)
        #print np.min(scattery), np.max(scattery)

        # First step is to set up the plot
        width = 20
        height = 5
        multiplier = 1.
        fig = plt.figure(figsize=(width * multiplier, height * multiplier))


        ax = fig.add_subplot(1, 3, 1) #, aspect='equal')
        ax.set_title('Residuals')
        ax.hist([loc[0]-loc[1] for loc in obs_sim_zone_all], bins=20, alpha=0.5)


        ax = fig.add_subplot(1, 3, 2) #, aspect=0.9)
        ax.set_title('Sim vs Obs (%d points)' %(len(scatterx)))

        comp_zone_plots={}        
        #colours = ['b', 'c', 'sienna', 'm', 'r', 'green', 'fuchsia']
        colours = ['r', 'orangered', 'y', 'green', 'teal', 'blue', 'fuchsia']
        for i in range(1,8):
            scatterx2 = [loc[0] for loc in obs_sim_zone_all if loc[2]==float(i)]        
            scattery2 = [loc[1] for loc in obs_sim_zone_all if loc[2]==float(i)]        
            #print len(scatterx2), colours[i-1]
            comp_zone_plots[i] = ax.scatter(scatterx2, scattery2, edgecolors=colours[i-1], facecolors='none', alpha=0.5)
        
        plt.legend((comp_zone_plots[1], comp_zone_plots[2], comp_zone_plots[3], 
                    comp_zone_plots[4], comp_zone_plots[5], comp_zone_plots[6],
                    comp_zone_plots[7]), 
                   ('qa', 'utb', 'utqa', 'utam', 'utaf', 'lta', 'bse'),
                   scatterpoints=1,
                   loc='upper left',
                   ncol=4,
                   fontsize=11)
        
        plt.xlabel('Observed')
        plt.ylabel('Simulated', labelpad=10)

        sum1 = 0.
        sum2 = 0.
        
        if len(scatterx) != 0:
        
            mean = np.mean(scatterx)
            for i in range(len(scatterx)):
                num1 = (scatterx[i] - scattery[i])          
                num2 = (scatterx[i] - mean)
                sum1 += num1 ** np.float64(2.)
                sum2 += num2 ** np.float64(2.)
            
            ME = 1 - sum1 / sum2
            
            ax.text(150, 75, 'Model Efficiency = %4.2f' %(ME))        
    
            # for PBIAS
            def pbias(simulated, observed):
                return np.sum(simulated-observed)*100/np.sum(observed)
    
            ax.text(150, 40, 'PBIAS = %4.2f%%' %(pbias(scattery, scatterx)))        
            
            # For rmse
            def rmse(simulated, observed):
                return np.sqrt(((simulated - observed) ** 2).mean())
    
            ax.text(150, 20, 'RMSE = %4.2f' %(rmse(scattery, scatterx)))        
        
        ax.plot(ax.get_ylim(), ax.get_ylim())        

        ax = fig.add_subplot(1, 3, 3) #, aspect=0.9)
        ax.set_title('Residuals in space')

        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        modelmap.plot_ibound()

        x = np.array([h[3] for h in obs_sim_zone_all])        
        y = np.array([h[4] for h in obs_sim_zone_all])        
        zone = [h[2] for h in obs_sim_zone_all]
        residuals = [h[0]-h[1] for h in obs_sim_zone_all]        
        residuals = np.absolute(residuals)

        #plt.scatter(x, y, c=residual, alpha=0.5)        

        #alphas = np.linspace(0.1, 1, 10)
        from matplotlib import colors
        import six
        colors_ = list(six.iteritems(colors.cnames))
        # Add the single letter colors.
        for name, rgb in six.iteritems(colors.ColorConverter.colors):
            hex_ = colors.rgb2hex(rgb)
            colors_.append((name, hex_))

        hex_ = [color[1] for color in colors_]
        nams = [color[0] for color in colors_]
        # Get the rgb equivalent.
        rgb_all = [colors.hex2color(color) for color in hex_]

        rgb_ref = []
        for col in colours:
            for index, nam in enumerate(nams):
                if col == nam:
                    print 'yep'
                    rgb_ref += [rgb_all[index]]
        print rgb_ref

        rgba_colors = np.zeros((len(x),4))
        # for red the first column needs to be one
        for i in range(1,8):
            rgba_colors[:,0][zone == i] = rgb_ref[i-1][0]
            rgba_colors[:,1][zone == i] = rgb_ref[i-1][1]
            rgba_colors[:,2][zone == i] = rgb_ref[i-1][2]
        # the fourth column needs to be your alphas
        rgba_colors[:, 3] = residuals / np.max(residuals) #alphas
        
        plt.scatter(x, y, color=rgba_colors)
        
        #fig.subplots_adjust(left=0.01, right=0.95, bottom=0.05, top=0.95, wspace=0.1, hspace=0.12)

        plt.show()
        
        
    def viewHeadsByZone(self, nper='all'):

        # Create the headfile object
        headobj = self.importHeads()
        times = headobj.get_times()        
        if nper == 'all':
            head = headobj.get_alldata()
            head = np.mean(head, axis=0)            
            head_orig = head
            zoned = self.HeadsByZone(head)
            head = zoned
        else:
            head = headobj.get_data(totim=times[nper])
            head_orig = head
            zoned = self.HeadsByZone(head)
            head = zoned

        if nper == 'all':
            scatterx = []
            scattery = []
            obs_sim_zone_all = []
            for i in range(len(times)):
                self.CompareObserved('head', head_orig, nper=i)                
                scatterx += [h[0] for h in self.obs_sim_zone]        
                scattery += [h[1] for h in self.obs_sim_zone]  
                obs_sim_zone_all += self.obs_sim_zone
            self.obs_sim_zone = obs_sim_zone_all
        else:
            self.CompareObserved('head', head_orig, nper=nper)                
            scatterx = [h[0] for h in self.obs_sim_zone]        
            scattery = [h[1] for h in self.obs_sim_zone]        

        # First step is to set up the plot
        width = 20
        height = 10
        multiplier = 1.
        fig = plt.figure(figsize=(width * multiplier, height * multiplier))

        vmin = 0
        vmax = 200

        ax = fig.add_subplot(2, 4, 1, aspect='equal')

        ax.set_title('ibound and bc')
        # Next we create an instance of the ModelMap class
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        modelmap.plot_ibound()
        #linecollection = modelmap.plot_grid()
        
        #ax = fig.add_subplot(1, 3, 2, aspect='equal')
        #ax.set_title('Riv BC')
        #modelmap.plot_bc(plotAll=True)        
        #modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        
        modelmap.plot_bc('RIV')
        modelmap.plot_bc('WEL')
        modelmap.plot_bc('GHB')
        #modelmap.plot_bc('DRN')
        ax.axes.xaxis.set_ticklabels([])

#        ax.set_title('qa')
#        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
#        max_head = np.amax(head)         
#        min_head = np.amin(head)
#        #print max_head
#        #print min_head
#        
#        array = modelmap.plot_array(head[0], masked_values=[-999.98999023, max_head, min_head], alpha=0.5, vmin=vmin, vmax=vmax)        
#        ax.xaxis.set_ticklabels([])
        #cbar_ax1 = fig.add_axes([0.19, 0.525, 0.01, 0.42])
        #fig.colorbar(array, cax=cbar_ax1)
        #linecollection = modelmap.plot_grid()
#        scatterx2 = [loc[1][0] for loc in self.obs_loc_val_zone if loc[2]==1.0]        
#        scattery2 = [loc[1][1] for loc in self.obs_loc_val_zone if loc[2]==1.0]        
#        ax.scatter(scatterx2, scattery2, c=[loc[0] for loc in self.obs_loc_val_zone if loc[2]==1.0], alpha=0.5, vmin=vmin, vmax=vmax)
#        ax.text(2.7e5, 6030000, 'Observations: %d' %(len(scatterx2)))        


        ax = fig.add_subplot(2, 4, 2, aspect='equal')
        ax.set_title('Coonambidgal')
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        max_head = np.amax(head)         
        min_head = np.amin(head)
        #print max_head
        #print min_head
        
        array = modelmap.plot_array(head[0], masked_values=[-999.98999023, max_head, min_head], alpha=0.5, vmin=vmin, vmax=vmax)        
#        ax.set_title('utb')
#        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
#        max_head = np.amax(head)         
#        min_head = np.amin(head)
#        #print max_head
#        #print min_head
#        
#        array = modelmap.plot_array(head[1], masked_values=[-999.98999023, max_head, min_head], alpha=0.5, vmin=vmin, vmax=vmax)        
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        cbar_ax2 = fig.add_axes([0.43, 0.525, 0.01, 0.42])
        fig.colorbar(array, cax=cbar_ax2)

        scatterx2 = [loc[3] for loc in self.obs_sim_zone if loc[2]==1.0]        
        scattery2 = [loc[4] for loc in self.obs_sim_zone if loc[2]==1.0]        
        ax.scatter(scatterx2, scattery2, c=[loc[0] for loc in self.obs_sim_zone if loc[2]==1.0], alpha=0.8, vmin=vmin, vmax=vmax)
#        scatterx2 = [loc[1][0] for loc in self.obs_sim_zone if loc[2] == 2.0]        
#        scattery2 = [loc[1][1] for loc in self.obs_sim_zone if loc[2] == 2.0]        
#        ax.scatter(scatterx2, scattery2, c=[loc[0] for loc in self.obs_sim_zone if loc[2] == 2.0], alpha=0.5, vmin=vmin, vmax=vmax)
        ax.text(2.7e5, 6030000, 'Observations: %d' %(len(scatterx2)))        


        ax = fig.add_subplot(2, 4, 3, aspect='equal')
        ax.set_title('Shepparton')
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        array = modelmap.plot_array(head[2], masked_values=[-999.98999023, max_head, min_head, np.nan], alpha=0.5, vmin=vmin, vmax=vmax)        
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        cbar_ax1 = fig.add_axes([0.67, 0.525, 0.01, 0.42])
        fig.colorbar(array, cax=cbar_ax1)

        scatterx2 = [loc[3] for loc in self.obs_sim_zone if loc[2] == 3.0]        
        scattery2 = [loc[4] for loc in self.obs_sim_zone if loc[2] == 3.0]        
        ax.scatter(scatterx2, scattery2, c=[loc[0] for loc in self.obs_sim_zone if loc[2] == 3.0], alpha=0.8, vmin=vmin, vmax=vmax)
        ax.text(2.7e5, 6030000, 'Observations: %d' %(len(scatterx2)))        


        ax = fig.add_subplot(2, 4, 4) #, aspect='equal')
        ax.set_title('Residuals')
        ax.hist([loc[0]-loc[1] for loc in self.obs_sim_zone], bins=20, alpha=0.5)
#        ax.set_title('utam')
#        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
#        array = modelmap.plot_array(head[3], masked_values=[-999.98999023, max_head, min_head, np.nan], alpha=0.5, vmin=vmin, vmax=vmax)        
#        ax.xaxis.set_ticklabels([])
#        ax.yaxis.set_ticklabels([])
#        cbar_ax5 = fig.add_axes([0.91, 0.525, 0.01, 0.42])
#        fig.colorbar(array, cax=cbar_ax5)
#
#        scatterx2 = [loc[1][0] for loc in self.obs_sim_zone if loc[2] == 4.0]        
#        scattery2 = [loc[1][1] for loc in self.obs_sim_zone if loc[2] == 4.0]        
#        ax.scatter(scatterx2, scattery2, c=[loc[0] for loc in self.obs_sim_zone if loc[2] == 4.0], alpha=0.5, vmin=vmin, vmax=vmax)
#        ax.text(2.7e5, 6030000, 'Observations: %d' %(len(scatterx2)))        
        
        ax = fig.add_subplot(2, 4, 5, aspect='equal')
        ax.set_title('Calivil')
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        array = modelmap.plot_array(head[4], masked_values=[-999.98999023, max_head, min_head, np.nan], alpha=0.5, vmin=vmin, vmax=vmax)        
        start, end = ax.get_xlim()
        start = start // 1000 * 1000 + 1000
        end = end // 1000 * 1000 - 1000
        ax.xaxis.set_ticks(np.arange(start, end, 20000.))
        ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
        cbar_ax3 = fig.add_axes([0.19, 0.055, 0.01, 0.42])
        fig.colorbar(array, cax=cbar_ax3)

        scatterx2 = [loc[3] for loc in self.obs_sim_zone if loc[2] == 5.0]        
        scattery2 = [loc[4] for loc in self.obs_sim_zone if loc[2] == 5.0]        
        ax.scatter(scatterx2, scattery2, c=[loc[0] for loc in self.obs_sim_zone if loc[2] == 5.0], alpha=0.8, vmin=vmin, vmax=vmax)
        ax.text(2.7e5, 6030000, 'Observations: %d' %(len(scatterx2)))        

        ax = fig.add_subplot(2, 4, 6, aspect='equal')
        ax.set_title('Renmark')
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        array = modelmap.plot_array(head[5], masked_values=[-999.98999023, max_head, min_head, np.nan], alpha=0.5, vmin=vmin, vmax=vmax)        
        ax.yaxis.set_ticklabels([])
        ax.xaxis.set_ticks(np.arange(start, end, 20000.))
        ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
        cbar_ax4 = fig.add_axes([0.43, 0.055, 0.01, 0.42])
        fig.colorbar(array, cax=cbar_ax4)

        scatterx2 = [loc[3] for loc in self.obs_sim_zone if loc[2] == 6.0]        
        scattery2 = [loc[4] for loc in self.obs_sim_zone if loc[2] == 6.0]        
        ax.scatter(scatterx2, scattery2, c=[loc[0] for loc in self.obs_sim_zone if loc[2] == 6.0], alpha=0.8, vmin=vmin, vmax=vmax)
        ax.text(2.7e5, 6030000, 'Observations: %d' %(len(scatterx2)))        

        ax = fig.add_subplot(2, 4, 7, aspect='equal')
        ax.set_title('Basement')
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        array = modelmap.plot_array(head[6], masked_values=[-999.98999023, max_head, min_head, np.nan], alpha=0.5, vmin=vmin, vmax=vmax)        
        ax.yaxis.set_ticklabels([])
        ax.xaxis.set_ticks(np.arange(start, end, 20000.))
        ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
        cbar_ax5 = fig.add_axes([0.67, 0.055, 0.01, 0.42])
        fig.colorbar(array, cax=cbar_ax5)

        scatterx2 = [loc[3] for loc in self.obs_sim_zone if loc[2] == 7.0]        
        scattery2 = [loc[4] for loc in self.obs_sim_zone if loc[2] == 7.0]        
        ax.scatter(scatterx2, scattery2, c=[loc[0] for loc in self.obs_sim_zone if loc[2] == 7.0], alpha=0.8, vmin=vmin, vmax=vmax)
        ax.text(2.7e5, 6030000, 'Observations: %d' %(len(scatterx2)))        

        ax = fig.add_subplot(2, 4, 8, aspect=0.9)
        ax.set_title('Sim vs Obs (%d points)' %(len(scatterx)))
        #ax.scatter(scatterx, scattery, edgecolor='none', c=[loc[2] for loc in self.obs_sim_zone], alpha=0.5, vmin=1, vmax=7)
        comp_zone_plots={}        
        colors = ['b', 'c', 'y', 'm', 'r', 'green', 'orange']
        for i in range(1,8):
            scatterx2 = [loc[0] for loc in self.obs_sim_zone if loc[2]==float(i)]        
            scattery2 = [loc[1] for loc in self.obs_sim_zone if loc[2]==float(i)]        
            comp_zone_plots[i] = ax.scatter(scatterx2, scattery2, edgecolors=colors[i-1], facecolors='none', alpha=0.5)
        
        plt.legend((comp_zone_plots[1], comp_zone_plots[2], comp_zone_plots[3], 
                    comp_zone_plots[4], comp_zone_plots[5], comp_zone_plots[6],
                    comp_zone_plots[7]), 
                   ('qa', 'utb', 'utqa', 'utam', 'utaf', 'lta', 'bse'),
                   scatterpoints=1,
                   loc='upper left',
                   ncol=4,
                   fontsize=11)
        
        plt.xlabel('Observed')
        plt.ylabel('Simulated', labelpad=10)

        scatterx = np.array(scatterx)        
        scattery = np.array(scattery)        
        sum1 = 0.
        sum2 = 0.
        
        if len(scatterx) != 0:
        
            mean = np.mean(scatterx)
            for i in range(len(scatterx)):
                num1 = (scatterx[i] - scattery[i])          
                num2 = (scatterx[i] - mean)
                sum1 += num1 ** np.float64(2.)
                sum2 += num2 ** np.float64(2.)
            
            ME = 1 - sum1 / sum2
            
            ax.text(150, 75, 'Model Efficiency = %4.2f' %(ME))        
    
            # for PBIAS
            def pbias(simulated, observed):
                return np.sum(simulated-observed)*100/np.sum(observed)
    
            ax.text(150, 40, 'PBIAS = %4.2f%%' %(pbias(scattery, scatterx)))        
            
            # For rmse
            def rmse(simulated, observed):
                return np.sqrt(((simulated - observed) ** 2).mean())
    
            ax.text(150, 20, 'RMSE = %4.2f' %(rmse(scattery, scatterx)))        
        
        ax.plot(ax.get_ylim(), ax.get_ylim())        
        #modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        #ffrf = modelmap.plot_array(head[6], masked_values=[-999.98999023, max_head, min_head], alpha=0.5)        
        
        #ax.yaxis.set_ticklabels([])
        #ax.xaxis.set_ticks(np.arange(start, end, 20000.))
        #ax.xaxis.get_major_formatter().set_powerlimits((0, 1))

        #cbar_ax5 = fig.add_axes([0.67, 0.055, 0.01, 0.42])
        #fig.colorbar(ffrf, cax=cbar_ax5)
        
        fig.subplots_adjust(left=0.01, right=0.95, bottom=0.05, top=0.95, wspace=0.1, hspace=0.12)

        plt.show()

    #End viewHeads        

    def viewHeadsByZone2(self, nper='all'):

        # Create the headfile object
        headobj = self.importHeads()
        times = headobj.get_times()        
        if nper == 'all':
            head = headobj.get_alldata()
            head = np.mean(head, axis=0)            
            head_orig = head
            zoned = self.HeadsByZone(head)
            head = zoned
        else:
            head = headobj.get_data(totim=times[nper])
            head_orig = head
            zoned = self.HeadsByZone(head)
            head = zoned

        if nper == 'all':
            scatterx = []
            scattery = []
            obs_sim_zone_all = []
            for i in range(len(times)):
                self.CompareObserved('head', head_orig, nper=i)                
                scatterx += [h[0] for h in self.obs_sim_zone]        
                scattery += [h[1] for h in self.obs_sim_zone]  
                obs_sim_zone_all += self.obs_sim_zone
            self.obs_sim_zone = obs_sim_zone_all
        else:
            self.CompareObserved('head', head_orig, nper=nper)                
            scatterx = [h[0] for h in self.obs_sim_zone]        
            scattery = [h[1] for h in self.obs_sim_zone]        

        # First step is to set up the plot
        width = 20
        height = 10
        multiplier = 1.
        fig = plt.figure(figsize=(width * multiplier, height * multiplier))

        vmin = 0
        vmax = 200

        ax = fig.add_subplot(2, 4, 1, aspect='equal')

        ax.set_title('ibound and bc')
        # Next we create an instance of the ModelMap class
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        modelmap.plot_ibound()
        #linecollection = modelmap.plot_grid()
        
        #ax = fig.add_subplot(1, 3, 2, aspect='equal')
        #ax.set_title('Riv BC')
        #modelmap.plot_bc(plotAll=True)        
        #modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        
        modelmap.plot_bc('RIV')
        modelmap.plot_bc('WEL')
        modelmap.plot_bc('GHB')
        #modelmap.plot_bc('DRN')
        #ax.axes.xaxis.set_ticklabels([])
        start, end = ax.get_xlim()
        start = start // 1000 * 1000 + 1000
        end = end // 1000 * 1000 - 1000
        ax.xaxis.get_major_formatter().set_powerlimits((0, 1))

#        ax.set_title('qa')
#        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
#        max_head = np.amax(head)         
#        min_head = np.amin(head)
#        #print max_head
#        #print min_head
#        
#        array = modelmap.plot_array(head[0], masked_values=[-999.98999023, max_head, min_head], alpha=0.5, vmin=vmin, vmax=vmax)        
#        ax.xaxis.set_ticklabels([])
        #cbar_ax1 = fig.add_axes([0.19, 0.525, 0.01, 0.42])
        #fig.colorbar(array, cax=cbar_ax1)
        #linecollection = modelmap.plot_grid()
#        scatterx2 = [loc[1][0] for loc in self.obs_loc_val_zone if loc[2]==1.0]        
#        scattery2 = [loc[1][1] for loc in self.obs_loc_val_zone if loc[2]==1.0]        
#        ax.scatter(scatterx2, scattery2, c=[loc[0] for loc in self.obs_loc_val_zone if loc[2]==1.0], alpha=0.5, vmin=vmin, vmax=vmax)
#        ax.text(2.7e5, 6030000, 'Observations: %d' %(len(scatterx2)))        


        ax = fig.add_subplot(2, 4, 2, aspect='equal')
        ax.set_title('Coonambidgal')
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        max_head = np.amax(head)         
        min_head = np.amin(head)
        #print max_head
        #print min_head
        
        array = modelmap.plot_array(head[0], masked_values=[-999.98999023, max_head, min_head], alpha=0.5, vmin=vmin, vmax=vmax)        
#        ax.set_title('utb')
#        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
#        max_head = np.amax(head)         
#        min_head = np.amin(head)
#        #print max_head
#        #print min_head
#        
#        array = modelmap.plot_array(head[1], masked_values=[-999.98999023, max_head, min_head], alpha=0.5, vmin=vmin, vmax=vmax)        
        #ax.xaxis.set_ticklabels([])
        #ax.yaxis.set_ticklabels([])
        ax.xaxis.set_ticks(np.arange(start, end, 20000.))
        ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
        cbar_ax2 = fig.add_axes([0.43, 0.525, 0.01, 0.42])
        fig.colorbar(array, cax=cbar_ax2)

        scatterx2 = [loc[3] for loc in self.obs_sim_zone if loc[2]==1.0]        
        scattery2 = [loc[4] for loc in self.obs_sim_zone if loc[2]==1.0]        
        ax.scatter(scatterx2, scattery2, c=[loc[0] for loc in self.obs_sim_zone if loc[2]==1.0], alpha=0.8, vmin=vmin, vmax=vmax)
#        scatterx2 = [loc[1][0] for loc in self.obs_sim_zone if loc[2] == 2.0]        
#        scattery2 = [loc[1][1] for loc in self.obs_sim_zone if loc[2] == 2.0]        
#        ax.scatter(scatterx2, scattery2, c=[loc[0] for loc in self.obs_sim_zone if loc[2] == 2.0], alpha=0.5, vmin=vmin, vmax=vmax)
        ax.text(2.7e5, 6030000, 'Observations: %d' %(len(scatterx2)))        


        ax = fig.add_subplot(2, 4, 3, aspect='equal')
        ax.set_title('Shepparton')
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        array = modelmap.plot_array(head[2], masked_values=[-999.98999023, max_head, min_head, np.nan], alpha=0.5, vmin=vmin, vmax=vmax)        
        #ax.xaxis.set_ticklabels([])
        #ax.yaxis.set_ticklabels([])
        ax.xaxis.set_ticks(np.arange(start, end, 20000.))
        ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
        cbar_ax1 = fig.add_axes([0.67, 0.525, 0.01, 0.42])
        fig.colorbar(array, cax=cbar_ax1)

        scatterx2 = [loc[3] for loc in self.obs_sim_zone if loc[2] == 3.0]        
        scattery2 = [loc[4] for loc in self.obs_sim_zone if loc[2] == 3.0]        
        ax.scatter(scatterx2, scattery2, c=[loc[0] for loc in self.obs_sim_zone if loc[2] == 3.0], alpha=0.8, vmin=vmin, vmax=vmax)
        ax.text(2.7e5, 6030000, 'Observations: %d' %(len(scatterx2)))        


        ax = fig.add_subplot(2, 4, 4) #, aspect='equal')
        ax.set_title('Residuals')
        ax.hist([loc[0]-loc[1] for loc in self.obs_sim_zone], bins=20, alpha=0.5)
#        ax.set_title('utam')
#        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
#        array = modelmap.plot_array(head[3], masked_values=[-999.98999023, max_head, min_head, np.nan], alpha=0.5, vmin=vmin, vmax=vmax)        
#        ax.xaxis.set_ticklabels([])
#        ax.yaxis.set_ticklabels([])
#        cbar_ax5 = fig.add_axes([0.91, 0.525, 0.01, 0.42])
#        fig.colorbar(array, cax=cbar_ax5)
#
#        scatterx2 = [loc[1][0] for loc in self.obs_sim_zone if loc[2] == 4.0]        
#        scattery2 = [loc[1][1] for loc in self.obs_sim_zone if loc[2] == 4.0]        
#        ax.scatter(scatterx2, scattery2, c=[loc[0] for loc in self.obs_sim_zone if loc[2] == 4.0], alpha=0.5, vmin=vmin, vmax=vmax)
#        ax.text(2.7e5, 6030000, 'Observations: %d' %(len(scatterx2)))        
        
        ax = fig.add_subplot(2, 4, 5, aspect='equal')
        ax.set_title('Calivil')
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        array = modelmap.plot_array(head[4], masked_values=[-999.98999023, max_head, min_head, np.nan], alpha=0.5, vmin=vmin, vmax=vmax)        
        ax.xaxis.set_ticks(np.arange(start, end, 20000.))
        ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
        cbar_ax3 = fig.add_axes([0.19, 0.055, 0.01, 0.42])
        fig.colorbar(array, cax=cbar_ax3)

        scatterx2 = [loc[3] for loc in self.obs_sim_zone if loc[2] == 5.0]        
        scattery2 = [loc[4] for loc in self.obs_sim_zone if loc[2] == 5.0]        
        ax.scatter(scatterx2, scattery2, c=[loc[0] for loc in self.obs_sim_zone if loc[2] == 5.0], alpha=0.8, vmin=vmin, vmax=vmax)
        ax.text(2.7e5, 6030000, 'Observations: %d' %(len(scatterx2)))        

        ax = fig.add_subplot(2, 4, 6, aspect='equal')
        ax.set_title('Renmark')
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        array = modelmap.plot_array(head[5], masked_values=[-999.98999023, max_head, min_head, np.nan], alpha=0.5, vmin=vmin, vmax=vmax)        
        #ax.yaxis.set_ticklabels([])
        ax.xaxis.set_ticks(np.arange(start, end, 20000.))
        ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
        cbar_ax4 = fig.add_axes([0.43, 0.055, 0.01, 0.42])
        fig.colorbar(array, cax=cbar_ax4)

        scatterx2 = [loc[3] for loc in self.obs_sim_zone if loc[2] == 6.0]        
        scattery2 = [loc[4] for loc in self.obs_sim_zone if loc[2] == 6.0]        
        ax.scatter(scatterx2, scattery2, c=[loc[0] for loc in self.obs_sim_zone if loc[2] == 6.0], alpha=0.8, vmin=vmin, vmax=vmax)
        ax.text(2.7e5, 6030000, 'Observations: %d' %(len(scatterx2)))        

        ax = fig.add_subplot(2, 4, 7, aspect='equal')
        ax.set_title('Basement')
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        array = modelmap.plot_array(head[6], masked_values=[-999.98999023, max_head, min_head, np.nan], alpha=0.5, vmin=vmin, vmax=vmax)        
        #ax.yaxis.set_ticklabels([])
        ax.xaxis.set_ticks(np.arange(start, end, 20000.))
        ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
        cbar_ax5 = fig.add_axes([0.67, 0.055, 0.01, 0.42])
        fig.colorbar(array, cax=cbar_ax5)

        scatterx2 = [loc[3] for loc in self.obs_sim_zone if loc[2] == 7.0]        
        scattery2 = [loc[4] for loc in self.obs_sim_zone if loc[2] == 7.0]        
        ax.scatter(scatterx2, scattery2, c=[loc[0] for loc in self.obs_sim_zone if loc[2] == 7.0], alpha=0.8, vmin=vmin, vmax=vmax)
        ax.text(2.7e5, 6030000, 'Observations: %d' %(len(scatterx2)))        

        ax = fig.add_subplot(2, 4, 8, aspect=0.9)
        ax.set_title('Sim vs Obs (%d points)' %(len(scatterx)))
        #ax.scatter(scatterx, scattery, edgecolor='none', c=[loc[2] for loc in self.obs_sim_zone], alpha=0.5, vmin=1, vmax=7)
        comp_zone_plots={}        
        colors = ['b', 'c', 'y', 'm', 'r', 'green', 'orange']
        for i in range(1,8):
            scatterx2 = [loc[0] for loc in self.obs_sim_zone if loc[2]==float(i)]        
            scattery2 = [loc[1] for loc in self.obs_sim_zone if loc[2]==float(i)]        
            comp_zone_plots[i] = ax.scatter(scatterx2, scattery2, edgecolors=colors[i-1], facecolors='none', alpha=0.5)
        
        plt.legend((comp_zone_plots[1], comp_zone_plots[2], comp_zone_plots[3], 
                    comp_zone_plots[4], comp_zone_plots[5], comp_zone_plots[6],
                    comp_zone_plots[7]), 
                   ('qa', 'utb', 'utqa', 'utam', 'utaf', 'lta', 'bse'),
                   scatterpoints=1,
                   loc='upper left',
                   ncol=4,
                   fontsize=11)
        
        plt.xlabel('Observed')
        plt.ylabel('Simulated', labelpad=10)

        scatterx = np.array(scatterx)        
        scattery = np.array(scattery)        
        sum1 = 0.
        sum2 = 0.
        
        if len(scatterx) != 0:
        
            mean = np.mean(scatterx)
            for i in range(len(scatterx)):
                num1 = (scatterx[i] - scattery[i])          
                num2 = (scatterx[i] - mean)
                sum1 += num1 ** np.float64(2.)
                sum2 += num2 ** np.float64(2.)
            
            ME = 1 - sum1 / sum2
            
            ax.text(150, 75, 'Model Efficiency = %4.2f' %(ME))        
    
            # for PBIAS
            def pbias(simulated, observed):
                return np.sum(simulated-observed)*100/np.sum(observed)
    
            ax.text(150, 40, 'PBIAS = %4.2f%%' %(pbias(scattery, scatterx)))        
            
            # For rmse
            def rmse(simulated, observed):
                return np.sqrt(((simulated - observed) ** 2).mean())
    
            ax.text(150, 20, 'RMSE = %4.2f' %(rmse(scattery, scatterx)))        
        
        ax.plot(ax.get_ylim(), ax.get_ylim())        
        #modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        #ffrf = modelmap.plot_array(head[6], masked_values=[-999.98999023, max_head, min_head], alpha=0.5)        
        
        #ax.yaxis.set_ticklabels([])
        #ax.xaxis.set_ticks(np.arange(start, end, 20000.))
        #ax.xaxis.get_major_formatter().set_powerlimits((0, 1))

        #cbar_ax5 = fig.add_axes([0.67, 0.055, 0.01, 0.42])
        #fig.colorbar(ffrf, cax=cbar_ax5)
        
        fig.subplots_adjust(left=0.01, right=0.95, bottom=0.05, top=0.95, wspace=0.1, hspace=0.12)

        plt.show()

    #End viewHeads    

    def viewHeads(self):

        # Create the headfile object
        headobj = self.importHeads()
        times = headobj.get_times()        
        head = headobj.get_data(totim=times[-1])

        #head_zoned = HeadsByZone(head)

        # First step is to set up the plot
        fig = plt.figure(figsize=(20, 10))
        ax = fig.add_subplot(2, 4, 1, aspect='equal')
        ax.set_title('ibound and bc')
        # Next we create an instance of the ModelMap class
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        modelmap.plot_ibound()
        #linecollection = modelmap.plot_grid()
        
        #modelmap.plot_bc('WEL')
        modelmap.plot_bc('RIV')
        ax.axes.xaxis.set_ticklabels([])

        #linecollection = modelmap.plot_grid()

        vmin = 0
        vmax = 200

        ax = fig.add_subplot(2, 4, 2, aspect='equal')
        ax.set_title('Heads layer 1')
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        max_head = np.amax(head)         
        min_head = np.amin(head)
        #print max_head
        #print min_head
        
        array = modelmap.plot_array(head[0], masked_values=[-999.98999023, max_head, min_head], alpha=0.5, vmin=vmin, vmax=vmax)        
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        cbar_ax2 = fig.add_axes([0.43, 0.525, 0.01, 0.42])
        fig.colorbar(array, cax=cbar_ax2)

        ax = fig.add_subplot(2, 4, 3, aspect='equal')
        ax.set_title('Heads layer 2')
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        array = modelmap.plot_array(head[1], masked_values=[-999.98999023, max_head, min_head], alpha=0.5, vmin=vmin, vmax=vmax)        
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        cbar_ax1 = fig.add_axes([0.67, 0.525, 0.01, 0.42])
        fig.colorbar(array, cax=cbar_ax1)

        ax = fig.add_subplot(2, 4, 4, aspect='equal')
        ax.set_title('Heads layer 3')
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        array = modelmap.plot_array(head[2], masked_values=[-999.98999023, max_head, min_head], alpha=0.5, vmin=vmin, vmax=vmax)        
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        cbar_ax5 = fig.add_axes([0.91, 0.525, 0.01, 0.42])
        fig.colorbar(array, cax=cbar_ax5)
        
        ax = fig.add_subplot(2, 4, 5, aspect='equal')
        ax.set_title('Heads layer 4')
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        array = modelmap.plot_array(head[3], masked_values=[-999.98999023, max_head, min_head], alpha=0.5, vmin=vmin, vmax=vmax)        
        start, end = ax.get_xlim()
        start = start // 1000 * 1000 + 1000
        end = end // 1000 * 1000 - 1000
        ax.xaxis.set_ticks(np.arange(start, end, 20000.))
        ax.xaxis.get_major_formatter().set_powerlimits((0, 1))

        cbar_ax3 = fig.add_axes([0.19, 0.055, 0.01, 0.42])
        fig.colorbar(array, cax=cbar_ax3)

        ax = fig.add_subplot(2, 4, 6, aspect='equal')
        ax.set_title('Heads layer 5')
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        array = modelmap.plot_array(head[4], masked_values=[-999.98999023, max_head, min_head], alpha=0.5, vmin=vmin, vmax=vmax)        
        ax.yaxis.set_ticklabels([])
        ax.xaxis.set_ticks(np.arange(start, end, 20000.))
        ax.xaxis.get_major_formatter().set_powerlimits((0, 1))

        cbar_ax4 = fig.add_axes([0.43, 0.055, 0.01, 0.42])
        fig.colorbar(array, cax=cbar_ax4)

        ax = fig.add_subplot(2, 4, 7, aspect='equal')
        ax.set_title('Heads layer 6')
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        array = modelmap.plot_array(head[5], masked_values=[-999.98999023, max_head, min_head], alpha=0.5, vmin=vmin, vmax=vmax)        
        ax.yaxis.set_ticklabels([])
        ax.xaxis.set_ticks(np.arange(start, end, 20000.))
        ax.xaxis.get_major_formatter().set_powerlimits((0, 1))

        cbar_ax5 = fig.add_axes([0.67, 0.055, 0.01, 0.42])
        fig.colorbar(array, cax=cbar_ax5)

    
        
        self.CompareObservedHead('head', head)                
        scatterx = [h[0] for h in self.comparison]        
        scattery = [h[1] for h in self.comparison]        

        ax = fig.add_subplot(2, 4, 8, aspect=0.9)
        ax.set_title('Sim vs Obs (%d points)' %(len(scatterx)))

                
        plt.scatter(scatterx, scattery, c=self.comp_zone)

        plt.xlabel('Observed')
        plt.ylabel('Simulated', labelpad=-20)
        
        scatterx = np.array(scatterx)        
        scattery = np.array(scattery)        
        sum1 = 0.
        sum2 = 0.
        mean = np.mean(scatterx)
        for i in range(len(scatterx)):
            num1 = (scatterx[i] - scattery[i])          
            num2 = (scatterx[i] - mean)
            sum1 += num1 ** np.float64(2.)
            sum2 += num2 ** np.float64(2.)
        
        ME = 1 - sum1 / sum2
        
        ax.text(150, -100, 'Model Efficiency = %4.2f' %(ME))        
        
        # For rmse
        def rmse(simulated, observed):
            return np.sqrt(((simulated - observed) ** 2).mean())

        ax.text(150, 0, 'RMSE = %4.2f' %(rmse(scattery, scatterx)))        
        
        ax.plot(ax.get_ylim(), ax.get_ylim())        
        #modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        #ffrf = modelmap.plot_array(head[6], masked_values=[-999.98999023, max_head, min_head], alpha=0.5)        
        
        #ax.yaxis.set_ticklabels([])
        #ax.xaxis.set_ticks(np.arange(start, end, 20000.))
        #ax.xaxis.get_major_formatter().set_powerlimits((0, 1))

        #cbar_ax5 = fig.add_axes([0.67, 0.055, 0.01, 0.42])
        #fig.colorbar(ffrf, cax=cbar_ax5)
        
        fig.subplots_adjust(left=0.01, right=0.95, bottom=0.05, top=0.95, wspace=0.1, hspace=0.12)

        plt.show()

    #End viewHeads        

    def viewHeads2(self):

        import flopy.utils.binaryfile as bf
        import pandas as pd
        #import matplotlib.ticker as plticker
        # Create the headfile object
        headobj = bf.HeadFile(self.data_folder + self.name+'.hds')
        cbbobj = bf.CellBudgetFile(self.data_folder + self.name+'.cbc')

        times = headobj.get_times()        
        head = headobj.get_data(totim=times[0])

        #print head
        #print dir(cbbobj)
        
        water_balance_components = cbbobj.textlist
        water_balance = {}        
        water_balance_summary = {}
        water_bal_summed_titles = []
        water_bal_summed_values = []
        for component in water_balance_components:
            component_stripped = component.lstrip().rstrip()
            if 'FLOW' in component_stripped: continue            
            water_balance[component_stripped] = cbbobj.get_data(text = component, full3D=True)[-1]            
            #print water_balance[component_stripped]            
            if np.any(water_balance[component_stripped][0]>0):
                water_balance_summary[component_stripped+'_pos'] = np.sum(water_balance[component_stripped][0][water_balance[component_stripped][0]>0])
            else:
                water_balance_summary[component_stripped+'_pos'] = 0.0
            if np.any(water_balance[component_stripped][0]<0):
                water_balance_summary[component_stripped+'_neg'] = np.sum(water_balance[component_stripped][0][water_balance[component_stripped][0]<0])
            else:
                water_balance_summary[component_stripped+'_neg'] = 0.0
            
            water_bal_summed_titles += component_stripped+'_pos', component_stripped+'_neg'
            water_bal_summed_values += water_balance_summary[component_stripped+'_pos'], water_balance_summary[component_stripped+'_neg']            
            #print component + ':'
            #print water_balance_summary[component_stripped+'_pos'], water_balance_summary[component_stripped+'_neg']
        
       #riv_flux = cbbobj.get_data(text = 'RIVER LEAKAGE', full3D=True)
        #wel_flux = cbbobj.get_data(text = 'WELLS', full3D=True)
        
        water_bal_summed_titles += ['Error']
        Error = sum(water_bal_summed_values)        
        water_bal_summed_values += [Error]       
        wat_bal_df = pd.DataFrame(water_bal_summed_values, water_bal_summed_titles) # pd.DataFrame.from_dict(water_balance_summary, orient='index')
        wat_bal_df.columns = ['Flux m^3/d']        
        wat_bal_df = wat_bal_df[wat_bal_df['Flux m^3/d'] != 0.]
        print wat_bal_df         

        #x = np.linspace(self.model_data.model_boundary[0], self.model_data.model_boundary[1], self.ncol)
        #y = np.linspace(self.model_data.model_boundary[2], self.model_data.model_boundary[3], self.nrow)
        #levels = np.arange(-200,200,10)
        #extent = (self.model_data.model_boundary[0] + self.delr/2., self.model_data.model_boundary[1] - self.delr/2., self.model_data.model_boundary[3] - self.delc/2., self.model_data.model_boundary[2] + self.delc/2.)
        #extent =  (self.model_data.model_boundary[0], self.model_data.model_boundary[1], self.model_data.model_boundary[3], self.model_data.model_boundary[2])    
        # First step is to set up the plot
        fig = plt.figure(figsize=(20, 10))
        ax = fig.add_subplot(2, 4, 1, aspect='equal')
        ax.set_title('ibound and bc')
        # Next we create an instance of the ModelMap class
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        modelmap.plot_ibound()
        #linecollection = modelmap.plot_grid()
        
        #ax = fig.add_subplot(1, 3, 2, aspect='equal')
        #ax.set_title('Riv BC')
        #modelmap.plot_bc(plotAll=True)        
        #modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        
        #modelmap.plot_bc('WEL')
        modelmap.plot_bc('RIV')
        try:        
            modelmap.plot_bc('GHB')
        except:
            pass
        ax.axes.xaxis.set_ticklabels([])

        #linecollection = modelmap.plot_grid()

        ax = fig.add_subplot(2, 4, 2, aspect='equal')
        ax.set_title('Heads')
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        max_head = np.amax(head)         
        min_head = np.amin(head)
        array = modelmap.plot_array(head, masked_values=[-999.98999023, max_head, min_head], alpha=0.5)        
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        cbar_ax2 = fig.add_axes([0.43, 0.525, 0.01, 0.42])
        fig.colorbar(array, cax=cbar_ax2)

        ax = fig.add_subplot(2, 4, 3, aspect='equal')
        ax.set_title('Riv exchange')
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        river_flux = modelmap.plot_array(water_balance['RIVER LEAKAGE'][0], alpha=0.5)        
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])

        cbar_ax1 = fig.add_axes([0.67, 0.525, 0.01, 0.42])
        fig.colorbar(river_flux, cax=cbar_ax1)
        
        # Setup params to get water balance aspect ratio looking nice
        aspect = float(12.5715/((wat_bal_df.max()[0]-wat_bal_df.min()[0])/float(wat_bal_df.shape[1])))        
        
        ax = fig.add_subplot(2, 4, 4, aspect=aspect)
        ax.set_title('Water Balance')
        wat_bal_df.plot(kind='bar', ax=plt.gca())        
        ax.grid(True)        
        gridlines = ax.get_xgridlines() #+ ax.get_ygridlines()
        for line in gridlines:
            line.set_linestyle('-')
    
        ax = fig.add_subplot(2, 4, 5, aspect='equal')
        ax.set_title('Rainfall recharge')
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        recharge = modelmap.plot_array(water_balance['RECHARGE'][0], masked_values=[0.], alpha=0.5)        
        #print water_balance['RECHARGE'][0]
        #print np.mean(water_balance['RECHARGE'][0])
        #print np.max(water_balance['RECHARGE'][0])
        #print np.min(water_balance['RECHARGE'][0])

        start, end = ax.get_xlim()
        start = start // 1000 * 1000 + 1000
        end = end // 1000 * 1000 - 1000
        ax.xaxis.set_ticks(np.arange(start, end, 20000.))
        ax.xaxis.get_major_formatter().set_powerlimits((0, 1))

        cbar_ax3 = fig.add_axes([0.19, 0.055, 0.01, 0.42])
        fig.colorbar(recharge, cax=cbar_ax3)

        ax = fig.add_subplot(2, 4, 6, aspect='equal')
        ax.set_title('Surface elevation')
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        elev = self.mf.dis.top.array
        #ibound_mask = np.ma.masked_where(self.model_data.model_mesh3D[1][0] == -1, self.model_data.model_mesh3D[1][0])
        elev = np.ma.masked_where(self.model_data.model_mesh3D[1][0] == 0, elev)
        elevation = modelmap.plot_array(elev, alpha=0.5)           
        ax.yaxis.set_ticklabels([])
        ax.xaxis.set_ticks(np.arange(start, end, 20000.))
        ax.xaxis.get_major_formatter().set_powerlimits((0, 1))

        cbar_ax4 = fig.add_axes([0.43, 0.055, 0.01, 0.42])
        fig.colorbar(elevation, cax=cbar_ax4)

        """
        ax = fig.add_subplot(2, 4, 6, aspect='equal')
        ax.set_title('GW Pumping')
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        pumping = modelmap.plot_array(water_balance['WELLS'][0], masked_values=[0.], alpha=0.5)        
        ax.yaxis.set_ticklabels([])
        ax.xaxis.set_ticks(np.arange(start, end, 20000.))
        ax.xaxis.get_major_formatter().set_powerlimits((0, 1))

        cbar_ax4 = fig.add_axes([0.43, 0.055, 0.01, 0.42])
        fig.colorbar(pumping, cax=cbar_ax4)
        """
        
        """        
        ax = fig.add_subplot(2, 4, 7, aspect='equal')
        ax.set_title('Storage change')
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        storage = modelmap.plot_array(water_balance['STORAGE'][0], masked_values=[0.], alpha=0.5)        
        ax.yaxis.set_ticklabels([])
        ax.xaxis.set_ticks(np.arange(start, end, 20000.))
        ax.xaxis.get_major_formatter().set_powerlimits((0, 1))

        cbar_ax5 = fig.add_axes([0.67, 0.055, 0.01, 0.42])
        fig.colorbar(storage, cax=cbar_ax5)
        """
        ax = fig.add_subplot(2, 4, 7, aspect='equal')
        ax.set_title('Pressure Head')
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        storage = modelmap.plot_array(head[0]-elev , masked_values=[0.], alpha=0.5, vmin=-50, vmax=50)        
        ax.yaxis.set_ticklabels([])
        ax.xaxis.set_ticks(np.arange(start, end, 20000.))
        ax.xaxis.get_major_formatter().set_powerlimits((0, 1))

        cbar_ax5 = fig.add_axes([0.67, 0.055, 0.01, 0.42])
        fig.colorbar(storage, cax=cbar_ax5)
        
        """
        ax = fig.add_subplot(2, 4, 8, aspect='equal')
        ax.set_title('Flow right face')
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        ffrf = modelmap.plot_array(water_balance['FLOW RIGHT FACE'][0], masked_values=[0.], alpha=0.5)        
        ax.yaxis.set_ticklabels([])
        ax.xaxis.set_ticks(np.arange(start, end, 20000.))
        ax.xaxis.set_major_formatter(plticker.FormatStrFormatter(fstring))

        cbar_ax5 = fig.add_axes([0.67, 0.055, 0.01, 0.42])
        fig.colorbar(ffrf, cax=cbar_ax5)
        """
        
        #plt.subplot(1, 1, 1, aspect='equal')
        #plt.title('stress period ' + str(iplot + 1))
        #plt.contourf(x, y, np.flipud(head[0, :, :]), cmap=plt.cm.rainbow, levels=levels, extent=extent)
        #plt.imshow(np.flipud(head[0, :, :]), cmap=plt.cm.rainbow)#, levels=levels, extent=extent)
        #plt.colorbar()
        
        #cb = plt.colorbar(array, shrink=0.5)        
        fig.subplots_adjust(left=0.01, right=0.95, bottom=0.05, top=0.95, wspace=0.1, hspace=0.12) #right=0.8

        #fig.subplots_adjust(right=0.8)


        #fig.tight_layout()
        plt.show()
        #plt.savefig('test_model.png')

        #self.mf.upw.hk.plot(masked_values=[0.], colorbar=True)
        #plt.show()

    #End viewHeads2        



#End ModflowModel()

class MT3DModel(object):
    
    def __init__(self):
        pass

    
    def createBTNpackage(self):
        #Add the BTN package to the model
        ibound = self.model_data.model_mesh3D[1]        
        ibound[ibound == -1] = 0
        self.btn = flopy.mt3d.Mt3dBtn(self.mt, icbund=ibound, ncomp=1, mcomp=1, 
                                      cinact=-9.9E1, thkmin=-1.0E-6, ifmtcn=5, 
                                      ifmtnp=0, ifmtrf=0, ifmtdp=0, nprs=0, 
                                      timprs=None, savucn=1, nprobs=0, 
                                      chkmas=1, nprmas=1, dt0=10000.0, ttsmax=100000.0)

    def createADVpackage(self):
        #Add the ADV package to the model
        self.adv = flopy.mt3d.Mt3dAdv(self.mt, mixelm=0, percel=1, 
                                      mxpart=250000, nadvfd=1, itrack=3, 
                                      wd=0.5, dceps=1.0E-4, nplane=0, npl=5, 
                                      nph=8, npmin=1, npmax=16, nlsink=0, 
                                      npsink=8, dchmoc=0.01)

    def createDSPpackage(self):
        #Add the DSP package to the model
        self.dfp = flopy.mt3d.Mt3dDsp(self.mt, multiDiff=True, al=10., trpt=0.1, trpv=0.1, dmcoef=0.0)

    def createRCTpackage(self):
        #Add the RCT package to the model
        self.rct = flopy.mt3d.Mt3dRct(self.mt, isothm=1, ireact=1, igetsc=0, rc1=np.log(2)/(5730*365))

    def createGCGpackage(self):
        #Add the GCG package to the model
        self.gcg = flopy.mt3d.Mt3dGcg(self.mt, mxiter=1000, iter1=100, isolve=1, 
                                      ncrs=0, accl=1, cclose=1.0E-4, iprgcg=0)
    
    def finaliseMT3Dmodel(self):
        self.mt.write_input()
    # end finaliseMT3Dmodel

    def buildMT3D(self):

        self.mt = flopy.mt3d.Mt3dms(modelname=self.name+'_transport', ftlfilename='mt3d_link.ftl', modflowmodel=self.mf, model_ws=self.data_folder, exe_name='MT3D-USGS_64.exe')

        self.createBTNpackage()
        self.createADVpackage()
        self.createDSPpackage()
        self.createRCTpackage()
        self.createGCGpackage()
        
        ssm_data = {}
        itype = flopy.mt3d.Mt3dSsm.itype_dict()
        for per in range(self.nper):
            ssm_data[per] = []
        ibound = self.model_data.model_mesh3D[1]        
        ibound[ibound == -1] = 0
        
        self.crch = {}
        for per in range(self.nper):
            self.crch[per] = []
        #end for
    
        for boundary in self.model_data.boundaries.bc:
            if self.model_data.boundaries.bc[boundary]['bc_type'] == 'recharge':
                for key in self.model_data.boundaries.bc[boundary]['bc_array'].keys():
                    self.crch[key] = np.ones_like(self.model_data.boundaries.bc[boundary]['bc_array'][key])               
                    self.crch[key] = self.crch[key] * 100.0
                    self.crch[key][ibound[0]==0] = 0.0
            #if self.model_data.boundaries.bc[boundary]['bc_type'] == 'river':
            #    self.createRIVpackage(self.model_data.boundaries.bc[boundary]['bc_array'])
                
            if self.model_data.boundaries.bc[boundary]['bc_type'] == 'wells':
                for key in self.model_data.boundaries.bc[boundary]['bc_array'].keys():
                    for well in self.model_data.boundaries.bc[boundary]['bc_array'][key]:
                        ssm_data[key].append((well[0], well[1], well[2], 100.0, itype['WEL']))
                        
            #if self.model_data.boundaries.bc[boundary]['bc_type'] == 'drain':
            #    self.model_data.boundaries.bc[boundary]['bc_array']

            if self.model_data.boundaries.bc[boundary]['bc_type'] == 'general head':
                self.model_data.boundaries.bc[boundary]['bc_array']
                for key in self.model_data.boundaries.bc[boundary]['bc_array'].keys():
                    for ghb in self.model_data.boundaries.bc[boundary]['bc_array'][key]:
                        ssm_data[key].append((ghb[0], ghb[1], ghb[2], 0.0, itype['GHB']))

        river_exists = False
        for boundary in self.model_data.boundaries.bc:
            if self.model_data.boundaries.bc[boundary]['bc_type'] == 'river':
                river_exists = True
                break
        
        if river_exists:
            river = {}
            river[0] = []
            for boundary in self.model_data.boundaries.bc:
                if self.model_data.boundaries.bc[boundary]['bc_type'] == 'river':
                    time_key = self.model_data.boundaries.bc[boundary]['bc_array'].keys()[0]
                    river[0] += self.model_data.boundaries.bc[boundary]['bc_array'][time_key]
                if self.model_data.boundaries.bc[boundary]['bc_type'] == 'channel':
                    time_key = self.model_data.boundaries.bc[boundary]['bc_array'].keys()[0]
                    river[0] += self.model_data.boundaries.bc[boundary]['bc_array'][time_key]
            
            for key in river.keys():
                for riv in river[key]:
                    ssm_data[key].append((riv[0], riv[1], riv[2], 100.0, itype['RIV']))

        self.ssm = flopy.mt3d.Mt3dSsm(self.mt, stress_period_data=ssm_data, crch=self.crch)

        self.finaliseMT3Dmodel()
    #End buildMODFLOW        

    def runMT3D(self):

        success, buff = self.mt.run_model()
        #if not success:
        #    raise Exception('MT3D did not terminate normally.')
        #End if
    #End runMT3D()    

    def importConcs(self):
        self.concobj = bf.UcnFile(self.data_folder + 'MT3D001.UCN')
        return self.concobj
   

    def viewConcsByZone(self, nper='all'):

        # Create the headfile object
        concobj = self.importConcs()
        times = concobj.get_times()        
        if nper == 'all':
            conc = concobj.get_alldata()
            conc = np.mean(conc, axis=0)            
            zoned = self.ConcsByZone(conc)
            conc = zoned
        else:
            conc = concobj.get_data(totim=times[nper])
            zoned = self.ConcsByZone(conc)
            conc = zoned

#        if nper == 'all':
#            scatterx = []
#            scattery = []
#            obs_sim_zone_all = []
#            for i in range(len(times)):
#                self.CompareObserved('head', head_orig, nper=i)                
#                scatterx += [h[0] for h in self.obs_sim_zone]        
#                scattery += [h[1] for h in self.obs_sim_zone]  
#                obs_sim_zone_all += self.obs_sim_zone
#            self.obs_sim_zone = obs_sim_zone_all
#        else:
#            self.CompareObserved('head', head_orig, nper=nper)                
#            scatterx = [h[0] for h in self.obs_sim_zone]        
#            scattery = [h[1] for h in self.obs_sim_zone]        

        # First step is to set up the plot
        width = 20
        height = 10
        multiplier = 1.
        fig = plt.figure(figsize=(width * multiplier, height * multiplier))

        vmin = 0.0
        vmax = 100.0

        ax = fig.add_subplot(2, 4, 1, aspect='equal')

        ax.set_title('ibound and bc')
        # Next we create an instance of the ModelMap class
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        modelmap.plot_ibound()
        #linecollection = modelmap.plot_grid()
        
        #ax = fig.add_subplot(1, 3, 2, aspect='equal')
        #ax.set_title('Riv BC')
        #modelmap.plot_bc(plotAll=True)        
        #modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        
        modelmap.plot_bc('RIV')
        modelmap.plot_bc('WEL')
        modelmap.plot_bc('GHB')
        #modelmap.plot_bc('DRN')
        ax.axes.xaxis.set_ticklabels([])

#        ax.set_title('qa')
#        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
#        max_head = np.amax(head)         
#        min_head = np.amin(head)
#        #print max_head
#        #print min_head
#        
#        array = modelmap.plot_array(head[0], masked_values=[-999.98999023, max_head, min_head], alpha=0.5, vmin=vmin, vmax=vmax)        
#        ax.xaxis.set_ticklabels([])
        #cbar_ax1 = fig.add_axes([0.19, 0.525, 0.01, 0.42])
        #fig.colorbar(array, cax=cbar_ax1)
        #linecollection = modelmap.plot_grid()
#        scatterx2 = [loc[1][0] for loc in self.obs_loc_val_zone if loc[2]==1.0]        
#        scattery2 = [loc[1][1] for loc in self.obs_loc_val_zone if loc[2]==1.0]        
#        ax.scatter(scatterx2, scattery2, c=[loc[0] for loc in self.obs_loc_val_zone if loc[2]==1.0], alpha=0.5, vmin=vmin, vmax=vmax)
#        ax.text(2.7e5, 6030000, 'Observations: %d' %(len(scatterx2)))        


        ax = fig.add_subplot(2, 4, 2, aspect='equal')
        ax.set_title('Coonambidgal')
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        max_conc = -100.0#np.amax(conc)         
        min_conc = 1E6#np.amin(conc)
        
        array = modelmap.plot_array(conc[0], masked_values=[-999.98999023, max_conc, min_conc], alpha=0.5, vmin=vmin, vmax=vmax)        
#        ax.set_title('utb')
#        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
#        max_head = np.amax(head)         
#        min_head = np.amin(head)
#        #print max_head
#        #print min_head
#        
#        array = modelmap.plot_array(head[1], masked_values=[-999.98999023, max_head, min_head], alpha=0.5, vmin=vmin, vmax=vmax)        
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        cbar_ax2 = fig.add_axes([0.43, 0.525, 0.01, 0.42])
        fig.colorbar(array, cax=cbar_ax2)

#        scatterx2 = [loc[3] for loc in self.obs_sim_zone if loc[2]==1.0]        
#        scattery2 = [loc[4] for loc in self.obs_sim_zone if loc[2]==1.0]        
#        ax.scatter(scatterx2, scattery2, c=[loc[0] for loc in self.obs_sim_zone if loc[2]==1.0], alpha=0.8, vmin=vmin, vmax=vmax)
#        ax.text(2.7e5, 6030000, 'Observations: %d' %(len(scatterx2)))        


        ax = fig.add_subplot(2, 4, 3, aspect='equal')
        ax.set_title('Shepparton')
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        array = modelmap.plot_array(conc[2], masked_values=[-999.98999023, max_conc, min_conc, np.nan], alpha=0.5, vmin=vmin, vmax=vmax)        
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        cbar_ax1 = fig.add_axes([0.67, 0.525, 0.01, 0.42])
        fig.colorbar(array, cax=cbar_ax1)

#        scatterx2 = [loc[3] for loc in self.obs_sim_zone if loc[2] == 3.0]        
#        scattery2 = [loc[4] for loc in self.obs_sim_zone if loc[2] == 3.0]        
#        ax.scatter(scatterx2, scattery2, c=[loc[0] for loc in self.obs_sim_zone if loc[2] == 3.0], alpha=0.8, vmin=vmin, vmax=vmax)
#        ax.text(2.7e5, 6030000, 'Observations: %d' %(len(scatterx2)))        


#        ax = fig.add_subplot(2, 4, 4) #, aspect='equal')
#        ax.set_title('Residuals')
#        ax.hist([loc[0]-loc[1] for loc in self.obs_sim_zone], bins=20, alpha=0.5)

#        ax.set_title('utam')
#        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
#        array = modelmap.plot_array(conc[3], masked_values=[-999.98999023, max_conc, min_conc, np.nan], alpha=0.5, vmin=vmin, vmax=vmax)        
#        ax.xaxis.set_ticklabels([])
#        ax.yaxis.set_ticklabels([])
#        cbar_ax5 = fig.add_axes([0.91, 0.525, 0.01, 0.42])
#        fig.colorbar(array, cax=cbar_ax5)
#
#        scatterx2 = [loc[1][0] for loc in self.obs_sim_zone if loc[2] == 4.0]        
#        scattery2 = [loc[1][1] for loc in self.obs_sim_zone if loc[2] == 4.0]        
#        ax.scatter(scatterx2, scattery2, c=[loc[0] for loc in self.obs_sim_zone if loc[2] == 4.0], alpha=0.5, vmin=vmin, vmax=vmax)
#        ax.text(2.7e5, 6030000, 'Observations: %d' %(len(scatterx2)))        
        
        ax = fig.add_subplot(2, 4, 5, aspect='equal')
        ax.set_title('Calivil')
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        array = modelmap.plot_array(conc[4], masked_values=[-999.98999023, max_conc, min_conc, np.nan], alpha=0.5, vmin=vmin, vmax=vmax)        
        start, end = ax.get_xlim()
        start = start // 1000 * 1000 + 1000
        end = end // 1000 * 1000 - 1000
        ax.xaxis.set_ticks(np.arange(start, end, 20000.))
        ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
        cbar_ax3 = fig.add_axes([0.19, 0.055, 0.01, 0.42])
        fig.colorbar(array, cax=cbar_ax3)

#        scatterx2 = [loc[3] for loc in self.obs_sim_zone if loc[2] == 5.0]        
#        scattery2 = [loc[4] for loc in self.obs_sim_zone if loc[2] == 5.0]        
#        ax.scatter(scatterx2, scattery2, c=[loc[0] for loc in self.obs_sim_zone if loc[2] == 5.0], alpha=0.8, vmin=vmin, vmax=vmax)
#        ax.text(2.7e5, 6030000, 'Observations: %d' %(len(scatterx2)))        

        ax = fig.add_subplot(2, 4, 6, aspect='equal')
        ax.set_title('Renmark')
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        array = modelmap.plot_array(conc[5], masked_values=[-999.98999023, max_conc, min_conc, np.nan], alpha=0.5, vmin=vmin, vmax=vmax)        
        ax.yaxis.set_ticklabels([])
        ax.xaxis.set_ticks(np.arange(start, end, 20000.))
        ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
        cbar_ax4 = fig.add_axes([0.43, 0.055, 0.01, 0.42])
        fig.colorbar(array, cax=cbar_ax4)

#        scatterx2 = [loc[3] for loc in self.obs_sim_zone if loc[2] == 6.0]        
#        scattery2 = [loc[4] for loc in self.obs_sim_zone if loc[2] == 6.0]        
#        ax.scatter(scatterx2, scattery2, c=[loc[0] for loc in self.obs_sim_zone if loc[2] == 6.0], alpha=0.8, vmin=vmin, vmax=vmax)
#        ax.text(2.7e5, 6030000, 'Observations: %d' %(len(scatterx2)))        

        ax = fig.add_subplot(2, 4, 7, aspect='equal')
        ax.set_title('Basement')
        modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        array = modelmap.plot_array(conc[6], masked_values=[-999.98999023, max_conc, min_conc, np.nan], alpha=0.5, vmin=vmin, vmax=vmax)        
        ax.yaxis.set_ticklabels([])
        ax.xaxis.set_ticks(np.arange(start, end, 20000.))
        ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
        cbar_ax5 = fig.add_axes([0.67, 0.055, 0.01, 0.42])
        fig.colorbar(array, cax=cbar_ax5)

#        scatterx2 = [loc[3] for loc in self.obs_sim_zone if loc[2] == 7.0]        
#        scattery2 = [loc[4] for loc in self.obs_sim_zone if loc[2] == 7.0]        
#        ax.scatter(scatterx2, scattery2, c=[loc[0] for loc in self.obs_sim_zone if loc[2] == 7.0], alpha=0.8, vmin=vmin, vmax=vmax)
#        ax.text(2.7e5, 6030000, 'Observations: %d' %(len(scatterx2)))        

#        ax = fig.add_subplot(2, 4, 8, aspect=0.9)
#        ax.set_title('Sim vs Obs (%d points)' %(len(scatterx)))
#        #ax.scatter(scatterx, scattery, edgecolor='none', c=[loc[2] for loc in self.obs_sim_zone], alpha=0.5, vmin=1, vmax=7)
#        comp_zone_plots={}        
#        colors = ['b', 'c', 'y', 'm', 'r', 'green', 'orange']
#        for i in range(1,8):
#            scatterx2 = [loc[0] for loc in self.obs_sim_zone if loc[2]==float(i)]        
#            scattery2 = [loc[1] for loc in self.obs_sim_zone if loc[2]==float(i)]        
#            comp_zone_plots[i] = ax.scatter(scatterx2, scattery2, edgecolors=colors[i-1], facecolors='none', alpha=0.5)
#        
#        plt.legend((comp_zone_plots[1], comp_zone_plots[2], comp_zone_plots[3], 
#                    comp_zone_plots[4], comp_zone_plots[5], comp_zone_plots[6],
#                    comp_zone_plots[7]), 
#                   ('qa', 'utb', 'utqa', 'utam', 'utaf', 'lta', 'bse'),
#                   scatterpoints=1,
#                   loc='upper left',
#                   ncol=4,
#                   fontsize=11)
#        
#        plt.xlabel('Observed')
#        plt.ylabel('Simulated', labelpad=10)
#
#        scatterx = np.array(scatterx)        
#        scattery = np.array(scattery)        
#        sum1 = 0.
#        sum2 = 0.
#        
#        if len(scatterx) != 0:
#        
#            mean = np.mean(scatterx)
#            for i in range(len(scatterx)):
#                num1 = (scatterx[i] - scattery[i])          
#                num2 = (scatterx[i] - mean)
#                sum1 += num1 ** np.float64(2.)
#                sum2 += num2 ** np.float64(2.)
#            
#            ME = 1 - sum1 / sum2
#            
#            ax.text(150, 75, 'Model Efficiency = %4.2f' %(ME))        
#    
#            # for PBIAS
#            def pbias(simulated, observed):
#                return np.sum(simulated-observed)*100/np.sum(observed)
#    
#            ax.text(150, 40, 'PBIAS = %4.2f%%' %(pbias(scattery, scatterx)))        
#            
#            # For rmse
#            def rmse(simulated, observed):
#                return np.sqrt(((simulated - observed) ** 2).mean())
#    
#            ax.text(150, 20, 'RMSE = %4.2f' %(rmse(scattery, scatterx)))        
#        
#        ax.plot(ax.get_ylim(), ax.get_ylim())        
        #modelmap = flopy.plot.ModelMap(model=self.mf) #, sr=self.mf.dis.sr, dis=self.mf.dis)
        #ffrf = modelmap.plot_array(conc[6], masked_values=[-999.98999023, max_conc, min_conc], alpha=0.5)        
        
        #ax.yaxis.set_ticklabels([])
        #ax.xaxis.set_ticks(np.arange(start, end, 20000.))
        #ax.xaxis.get_major_formatter().set_powerlimits((0, 1))

        #cbar_ax5 = fig.add_axes([0.67, 0.055, 0.01, 0.42])
        #fig.colorbar(ffrf, cax=cbar_ax5)
        
        fig.subplots_adjust(left=0.01, right=0.95, bottom=0.05, top=0.95, wspace=0.1, hspace=0.12)

        plt.show()

#End viewConcsByZone      
    #End MT3DModel
        

if __name__ == '__main__':
    pass

#End main
