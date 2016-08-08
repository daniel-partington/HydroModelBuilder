import multiprocessing
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
        self.nper = 1
        self.perlen = 1 #36500 
        self.nstp = 1 #10
        self.steady = False
        self.start_datetime = "1/1/1970"
        
        # Initial data:
        self.strt = self.model_data.initial_conditions.ic_data["Head"] #self.model_data.model_mesh3D[0][1:] + 20.         
        
        self.hk = self.model_data.model_mesh3D[1].astype(float)
        self.hk[self.hk > 0] = 10.0
        self.hk[self.hk == -1] = 1.
        self.vka = 1.0 
        self.sy = 0.01 
        self.ss = 1.0E-4        

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
        
        strt = self.strt #np.nan*np.empty((nlay, nrow, ncol), dtype=np.float32)

        self.bas = flopy.modflow.ModflowBas(self.mf, ibound=ibound, strt=strt)

    #End setupBASPackage()

    def setupNWTpackage(self):
        # Specify NWT settings
        self.nwt = flopy.modflow.ModflowNwt(self.mf, headtol=1.0E1, fluxtol=1.0E1, thickfact=1E-7, options='COMPLEX')

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

    def finaliseModel(self):
        # Write the MODFLOW model input files
        self.mf.write_input()        

    # end finaliseModel

    def buildMODFLOW(self):
        pass
    #End buildMODFLOW        

    def runMODFLOW(self):

        self.name += str(multiprocessing.current_process().name) # so each process has own files
        self.mf = flopy.modflow.Modflow(self.name, exe_name=self.executable, model_ws=self.data_folder, version='mfnwt')

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
                
            if self.model_data.boundaries.bc[boundary]['bc_type'] == 'wells':
                self.createWELpackage(self.model_data.boundaries.bc[boundary]['bc_array'])

        river = {}
        river[0] = []
        for boundary in self.model_data.boundaries.bc:
            if self.model_data.boundaries.bc[boundary]['bc_type'] == 'river':
                time_key = self.model_data.boundaries.bc[boundary]['bc_array'].keys()[0]
                river[0] += self.model_data.boundaries.bc[boundary]['bc_array'][time_key]
        
        self.createRIVpackage(river)
        
        self.finaliseModel()

        # Run the MODFLOW model
        success, buff = self.mf.run_model()
        if not success:
            raise Exception('MODFLOW did not terminate normally.')
        #End if

    #End runMODFLOW()

    def HeadsByZone(self, heads):
        
        self.model_data.model_mesh3D[1]

        heads_zoned = np.full(self.model_data.model_mesh3D[1].shape, np.Nan)
        
        for zone in range(np.max(self.model_data.model_mesh3D[1])):
            
            heads_zoned[zone] = heads[zone]

        return heads_zoned

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
                        times = headobj.get_times()
                        head = headobj.get_data(totim=times[0])
                elif obs_type == 'stage':
                    #stage = self.importStage()
                    pass
                #end if
                obs_df = self.model_data.observations.obs_group[obs_set]['time_series']
                obs_df = obs_df[obs_df['active'] == True]
                sim_map_dict = self.model_data.observations.obs_group[obs_set]['mapped_observations']        
                for observation in obs_df['name']:
                    sim_head = head[sim_map_dict[observation][0]][sim_map_dict[observation][1]][sim_map_dict[observation][2]]
                    f.write('%f\n' %sim_head)

    def CompareObservedHead(self, obs_set, head):
        
        self.comparison = []
        self.comp_zone = []
        self.obs_loc_val_lay = []
        obs_df = self.model_data.observations.obs_group[obs_set]['time_series']
        obs_df = obs_df[obs_df['active'] == True]
        sim_map_dict = self.model_data.observations.obs_group[obs_set]['mapped_observations']        
        for observation in obs_df['name']: #self.model_data.observations.obs_group[obs_set]['time_series']['name']:
            idx = obs_df[obs_df['name'] == observation].index.tolist()[0]            
            obs = obs_df.get_value(idx, 'value')
            sim = head[sim_map_dict[observation][0]][sim_map_dict[observation][1]][sim_map_dict[observation][2]]

            # The following commented block should no longer be required as these
            # conditions should no longer arise due to prefiltering of observations.
#            if sim == np.float32(-999.99):
#                continue
#            if sim == np.float32(1E+30): 
#                continue

            self.comparison += [[obs, sim]]    
            self.comp_zone += self.model_data.model_mesh3D[1][sim_map_dict[observation][0]][sim_map_dict[observation][1]][sim_map_dict[observation][2]]/np.max(self.model_data.model_mesh3D[1])
            
            #x = self.model_data.observations.obs[obs_set]['locations']            
            #y = self.model_data.observations.obs[obs_set]['locations']            
            #self.obs_loc_val_lay += [[obs, x, y, sim_map_dict[observation][0]]]
            
    def importHeads(self):
        self.headobj = bf.HeadFile(self.data_folder + self.name+'.hds')
        return self.headobj

    def importCbb(self):
        self.cbbobj = bf.CellBudgetFile(self.data_folder + self.name+'.cbc')
        return self.cbbobj
      
    def viewHeads(self):

        # Create the headfile object
        headobj = self.importHeads()
        times = headobj.get_times()        
        head = headobj.get_data(totim=times[0])

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
        plt.ylabel('Simulated', labelpad=-10)
        
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
            water_balance[component_stripped] = cbbobj.get_data(text = component, full3D=True)            
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
        start, end = ax.get_xlim()
        start = start // 1000 * 1000 + 1000
        end = end // 1000 * 1000 - 1000
        ax.xaxis.set_ticks(np.arange(start, end, 20000.))
        ax.xaxis.get_major_formatter().set_powerlimits((0, 1))

        cbar_ax3 = fig.add_axes([0.19, 0.055, 0.01, 0.42])
        fig.colorbar(recharge, cax=cbar_ax3)

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

    def getRiverAquiferFlux(self, from_cell, to_cell):
        """
        Function to retrieve the flux between "from_cell" and "to_cell", works
        on the assumption that the stream cells are in order from top of stream
        to the bottom of the stream and includes fluxes for the bounding cells.
        """
        exchange = 0.0
        riv_flux = self.cbbobj.get_data(text = 'RIVER LEAKAGE', full3D=True)
        
        for cell in reach_cells:
            exchange += 1
        
        return exchange


#End Groundwater()

#Ordinarily, the script and the Class definition would be in separate files
if __name__ == '__main__':
    pass

#End main
