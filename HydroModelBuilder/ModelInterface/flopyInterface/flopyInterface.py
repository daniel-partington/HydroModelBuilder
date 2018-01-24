import datetime
import inspect
import os
import warnings
from types import MethodType

import flopy
import flopy.utils.binaryfile as bf
import numpy as np
import pandas as pd
from flopy.utils.sfroutputfile import SfrFile

import viz.flopy_viz as fviz  # Visualization extension methods
from Radon_EC_simple import Radon_EC_simple


class ModflowModel(object):

    def __init__(self, model_data, data_folder=None, **kwargs):
        """
        :param model_data: ModelManager instance, containing all the data for the model
        :param data_folder: str, path to data folder (MODFLOW model run files)
        """
        self.model_data = model_data
        self.name = model_data.name
        if data_folder is None:
            # the addition of os.path.sep is for legacy purposes
            self.data_folder = os.path.join(os.getcwd(), 'model_' + self.name) + os.path.sep
        else:
            if not os.path.exists(data_folder):
                # Creating temporary run folder
                # this may cause a race condition - if the directory is created/removed between check and creation
                # in parallel applications
                os.makedirs(data_folder)
            self.data_folder = data_folder + os.path.sep
        # End if

        self.executable = r".\MODFLOW-NWT_64.exe"

        self.nlay = self.model_data.model_mesh3D[0].shape[0] - 1
        self.nrow = self.model_data.model_mesh3D[0].shape[1]
        self.ncol = self.model_data.model_mesh3D[0].shape[2]
        self.delr = self.model_data.gridHeight
        self.delc = self.model_data.gridWidth
        self.top = self.model_data.model_mesh3D[0][0]
        self.botm = self.model_data.model_mesh3D[0][1:]
        self.xul = self.model_data.model_boundary[0]
        self.yul = self.model_data.model_boundary[3]

        # Need to document why these values are used
        self.headtol = 1E-3  # 1E-6
        self.fluxtol = 1.0E4

        if self.model_data.model_time.t['steady_state']:
            self.nper = 1
            self.perlen = 1
            self.perlen = 40000 * 365
            self.nstp = 1  # 0
            self.steady = True
        else:
            self.nper = self.model_data.model_time.t['steps']
            self.perlen = [
                x.total_seconds() / 86400.0 for x in self.model_data.model_time.t['intervals']]
            self.nstp = 1
            self.steady = False
            self.start_datetime = self.model_data.model_time.t['dateindex'][0]

        # Initial data:
        self.strt = self.model_data.initial_conditions.ic_data["Head"]

        self.hk = self.model_data.properties.properties['Kh']
        self.hk[self.hk == -1] = 1.0
        self.vka = self.model_data.properties.properties['Kv']
        self.sy = self.model_data.properties.properties['Sy']
        self.ss = self.model_data.properties.properties['SS']

        # Set all other kwargs as class attributes
        for key, value in kwargs.items():
            setattr(self, key, value)
        # End For

        # Add visualization methods to this instance
        viz_methods = inspect.getmembers(fviz, inspect.isfunction)
        for method_name, func in viz_methods:
            setattr(self, method_name, MethodType(func, self))
        # End for

    # End __init__()

    def createDiscretisation(self):
        # Create discretisation object
        self.dis = flopy.modflow.ModflowDis(self.mf,
                                            nlay=self.nlay,
                                            nrow=self.nrow,
                                            ncol=self.ncol,
                                            delr=self.delr,
                                            delc=self.delc,
                                            top=self.top,
                                            botm=self.botm,
                                            xul=self.xul,
                                            yul=self.yul,
                                            nper=self.nper,
                                            perlen=self.perlen,
                                            nstp=self.nstp,
                                            steady=self.steady)
        if self.verbose and self.check:
            self.dis.check()
        # End if
    # End createDiscretisation()

    def setupBASPackage(self, nlay, nrow, ncol):
        """
        Create and setup MODFLOW BAS package.

        :param nlay: int, number of layers
        :param nrow: int, number of rows
        :param ncol: int, number of columns
        """
        # Variables for the BAS package
        ibound = self.model_data.model_mesh3D[1]
        ibound[ibound == -1] = 0

        self.bas = flopy.modflow.ModflowBas(self.mf, ibound=ibound, strt=self.strt)
        if self.verbose and self.check:
            self.bas.check()
        # End if
    # End setupBASPackage()

    def setupNWTpackage(self, headtol, fluxtol):
        # Specify NWT settings
        self.nwt = flopy.modflow.ModflowNwt(self.mf,
                                            headtol=headtol,
                                            fluxtol=fluxtol,
                                            linmeth=2,
                                            iprnwt=1,
                                            ibotav=1,
                                            thickfact=1E-5,
                                            maxiterout=10000,
                                            options='COMPLEX')
    # End setupNWTpackage()

    def setupPCGpackage(self):
        self.pcg = flopy.modflow.ModflowPcg(self.mf)  # if using mf 2005
    # End setupPCGpachage()

    def setupUPWpackage(self):
        """
        This function is used to set up the upstream weighting package for the
        NWT solver in MODFLOW
        """
        # Add UPW package to the MODFLOW model to represent aquifers
        self.upw = flopy.modflow.ModflowUpw(self.mf,
                                            hk=self.hk,
                                            # layvka=self.vka_ani_flag
                                            vka=self.vka,
                                            sy=self.sy,
                                            ss=self.ss,
                                            hdry=-999.9,
                                            laywet=0,
                                            laytyp=1)

    # End setupUPWpackage()

    def createRIVpackage(self, lrcd=None):
        self.riv = flopy.modflow.ModflowRiv(self.mf, ipakcb=53, stress_period_data=lrcd)
        if self.verbose and self.check:
            self.riv.check()
        # End if
    # End createRIVpackage()

    def createSTRpackage(self, STRvariables=None):
        self.str = flopy.modflow.ModflowStr(self.mf, ipakcb=53, stress_period_data=STRvariables)
        if self.verbose and self.check:
            self.str.check()
    # End createSTRpackage()

    def createSFRpackage(self, reach_data, segment_data):
        if self.model_data.model_time.t['steady_state']:
            transroute = False
        else:
            transroute = True
        nstrm = reach_data.shape[0]
        nss = segment_data[0].shape[0]
        dataset_5 = {key: [nss, 0, 0] for key in segment_data.keys()}
        if self.verbose:
            print('SFR: Assuming units of model are m^3/d, change "const" if this is not true')
        # End if
        self.sfr = flopy.modflow.ModflowSfr2(self.mf,
                                             nstrm=nstrm,  # number of stream reaches, stream cells
                                             nss=nss,  # number of stream segments
                                             nsfrpar=0,  # DO NOT CHANGE
                                             nparseg=0,  # Using reach input so set this to zero
                                             const=86400,  # This is 1.0 for m^3/s and otherwise if units different
                                             dleak=0.0001,  # Tolerance level of stream depth used in computing leakage between each stream reach and active model cell
                                             ipakcb=53,  # Write leakage to cell by cell file
                                             istcb2=81,  # Write flow info to output file
                                             isfropt=1,
                                             # Next three only used if isfrop > 1
                                             nstrail=10,  # IGNORE FOR NOW
                                             isuzn=1,  # IGNORE FOR NOW
                                             nsfrsets=30,  # IGNORE FOR NOW
                                             # This can only be 1 at the moment and specifies the use of the kinematic
                                             # wave (KW) equation for solving stream flow
                                             irtflg=1,
                                             numtim=1,  # Number of timesteps to use within each MF-NWT timestep
                                             weight=0.75,  # Time weighting factor to calculate the change in channel storage
                                             flwtol=0.00003,  # flow tolerance convergence criterion for KW equation
                                             reach_data=reach_data,  # "Data Set 2", rec array with as many entries as stream reaches
                                             segment_data=segment_data,  # "Data Set 4b", dict with recarray of segment data for each stress period
                                             channel_geometry_data=None,  # "Data Set 6d", Can ignore for now, but necessary if icalc is > 1
                                             channel_flow_data=None,  # "Data Set 6e"
                                             dataset_5=dataset_5,  # This is not required unless one wants to mod printing across stress periods
                                             reachinput=True,
                                             transroute=transroute,
                                             tabfiles=False,
                                             tabfiles_dict=None,
                                             extension='sfr',  # Leave as is which is the default
                                             unit_number=None,  # Leave as None and let flopy define unit_number
                                             filenames=None  # Leave as None and let flopy define output names
                                             )

        if self.verbose and self.check:
            self.sfr.check()
        # End if
    # End createSFRpackage()

    def createGagepackage(self, gages, files=None):
        # gages should contain seg, rch, unit, outtype [set = 9]
        if files is None:
            files = ['Gage.gage']
            files += ['Gage.gag{}'.format(x) for x in range(len(gages))]
        # End if

        self.gage = flopy.modflow.ModflowGage(self.mf, numgage=len(gages), gage_data=gages, filenames=None)
        if self.verbose and self.check:
            self.gage.check()
        # End if
    # End createGagepackage()

    def createDRNpackage(self, lrcsc=None):
        self.drn = flopy.modflow.ModflowDrn(self.mf, ipakcb=53, stress_period_data=lrcsc)
        if self.verbose and self.check:
            self.drn.check()
    # End createDRNpackage()

    def createGHBpackage(self, lrcsc=None):
        self.ghb = flopy.modflow.ModflowGhb(self.mf, ipakcb=53, stress_period_data=lrcsc)
        if self.verbose and self.check:
            self.ghb.check()
    # end createGHBpackage

    def createRCHpackage(self, rchrate=None):
        """
        Add RCH package to the MODFLOW model to represent recharge
        """
        self.rch = flopy.modflow.ModflowRch(self.mf, ipakcb=53, rech=rchrate, nrchop=3)
        if self.verbose and self.check:
            self.rch.check()
    # end createRCHpackage

    def createWELpackage(self, lrcq=None):
        """
        Add WEL package to the MODFLOW model to represent pumping wells
        Expects a dictionary with an array of well location lay row col and flux
        at each stress period, e.g.

        >>> lrcq = {}
        >>> lrcq[0] = [[0, 7, 7, -100.]] # layer, row, column, flux
        """
        self.wel = flopy.modflow.ModflowWel(self.mf, ipakcb=53, stress_period_data=lrcq)
        if self.verbose and self.check:
            self.wel.check()
    # end createWElpackage

    def createOCpackage(self):
        """
        Add OC package to the MODFLOW model
        """
        spd_opts = ['save head', 'print budget', 'save budget']
        if self.steady:
            spd = {(0, 0): spd_opts[:],
                   (0, 1): spd_opts[:]}
        else:
            spd = {}
            for i in range(self.nper):
                spd[(i, 0)] = spd_opts[:]
                spd[(i, 1)] = spd_opts[:]
            # End for
        # End if

        self.oc = flopy.modflow.ModflowOc(self.mf, stress_period_data=spd, cboufm='(20i5)')
    # End createOCpackage()

    def createLMTpackage(self):
        """
        Add LMT package to the MODFLOW model to allow linking with MT3DMS
        """
        self.lmt = flopy.modflow.ModflowLmt(self.mf,
                                            output_file_header='extended',
                                            output_file_format='formatted',
                                            package_flows=['sfr'])
    # End createLMTpackage()

    def finaliseModel(self):
        # Write the MODFLOW model input files
        warnings.warn("Deprecated method called. Use `finalize_model()` instead", DeprecationWarning)
        self.finalize_model()
    # end finaliseModel

    def finalize_model(self):
        self.mf.write_input()
    # End finalize_model()

    def add_bc_to_target(self, bc_array, target):
        for key in bc_array:
            try:
                target[key] += bc_array[key]
            except KeyError:
                target[key] = bc_array[key]
            # End try
        # End for

        return target
    # End add_bc_to_river()

    def buildMODFLOW(self, transport=False, write=True, verbose=True, check=False):
        """
        Build MODFLOW model.
        """
        self.mf = flopy.modflow.Modflow(self.name, exe_name=self.executable,
                                        model_ws=self.data_folder, version='mfnwt')
        self.verbose = verbose
        self.check = check

        self.createDiscretisation()
        self.setupBASPackage(self.nlay, self.nrow, self.ncol)
        self.setupNWTpackage(self.headtol, self.fluxtol)
        self.setupUPWpackage()
        self.createOCpackage()

        river_exists = False
        wells_exist = False
        river = {}
        wel = {}
        bc = self.model_data.boundaries.bc
        for boundary in bc:
            bc_boundary = bc[boundary]
            bc_type = bc_boundary['bc_type']
            bc_array = bc_boundary['bc_array']

            if bc_type == 'recharge':
                self.createRCHpackage(rchrate=bc_array)

            if bc_type == 'drain':
                self.createDRNpackage(bc_array)

            if bc_type == 'general head':
                self.createGHBpackage(bc_array)

            if (bc_type == 'river') or (bc_type == 'channel'):
                river_exists = True
                river = self.add_bc_to_target(bc_array, river)
            # End if

            if (bc_type == 'river_flow'):
                self.createSFRpackage(bc_array[0], bc_array[1])
            # End if

            if bc_type == 'wells':
                wells_exist = True
                wel = self.add_bc_to_target(bc_array, wel)
            # End if
        # End for

        if river_exists:
            self.createRIVpackage(river)

        if wells_exist:
            self.createWELpackage(wel)

        # IF transport then
        if transport:
            self.createLMTpackage()

        if write:
            self.finalize_model()

        if self.verbose and self.check:
            self.checkMODFLOW()
        # End if

    # End buildMODFLOW()

    def checkMODFLOW(self):
        self.mf.check()
    # End checkMODFLOW

    def runMODFLOW(self, silent=True):
        '''
        Function to run the modflow model with the optional argument for 
        suppressing any output to screen.
        
        Call to checkConvergence is to capture model run failures that are not
        picked up in the flopy run_model function.
        '''

        success, buff = self.mf.run_model(silent=silent)
<<<<<<< Updated upstream

        return self.checkConvergence(fail=not success)

=======
        return success
>>>>>>> Stashed changes
    # End runMODFLOW()

    def checkConvergence(self, path=None, name=None, fail=False):
        converge_fail_options = ["****FAILED TO MEET SOLVER CONVERGENCE CRITERIA IN TIME STEP",  # Clear statement of model fail in list file
                                 " PERCENT DISCREPANCY =         200.00",  # Convergence but extreme discrepancy in results
                                 " NaN "  # Something big went wrong but somehow convergence was reached?
                                 ]
        path = path if path else self.data_folder
        name = name if name else self.name
        with open(os.path.join(path, name + '.list'), 'r') as f:
            list_file = f.read()

        # TODO consolidate these file writing processes, see `Utilities.text_writer`
        for converge_fail in converge_fail_options:
            if converge_fail in list_file:
                print("*** Convergence failure ***")
                now = datetime.datetime.now().strftime("%I%M%p%B%d%Y")
                with open(os.path.join(self.data_folder, "converge_fail_%s.txt" % now), 'w') as f:
                    f.write("Model did not converge, @ %s" %
                               datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y"))
                    f.write("Error: \n {}".format(converge_fail))
                return False
            # End if

            if fail:
                print("*** Model run failure ***")
                now = datetime.datetime.now().strftime("%I%M%p%B%d%Y")
                with open(os.path.join(self.data_folder, "converge_fail_%s.txt" % now), 'w') as f:
                    f.write("Model did not run, @ %s" %
                               datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y"))
                    f.write("Error: \n {}".format('MODFLOW did not terminate normally'))
                return False
            # End if
        # End for

        return True
    # End checkConvergence()

    def Calculate_Rn_from_SFR_with_simple_model(self, df, Ini_cond, Rn_decay=0.181):
        '''
        Use a simple model to calculate Radon concentrations in the stream
        based on outputs from the SFR package and using some other data that
        is required as arguments to this function.

        In particular the df is a pandas dataframe that should contain the output
        from the sfr package which can be imported via the sfroutputfile util
        from flopy. The other parameters required are:

        * Ini_Cond = Ini_cond # 3 item list containing Initial flow, radon and ec concentrations
        * Rn_decay = Rn_decay # constant for Radon decay

        # Dataframe variables
        * df['HZ_poro'] # Hyporheic zone porosity
        * df['HZ_Prod_Rate'] # Hyporheic zone radon production rate
        * df['HZ_RTime'] # Hyporheic zone residence time
        * df['R_depth_HZ'] # Depth of the hyporheic zone
        * df['GTV'] # Gas transfer velocity
        * df['GW_Rn_conc'] # Groundwater radon concentration
        * df['GW_EC'] # Groundwater EC
        * df['Trib_EC'] # EC of the inflowing tributary water if present
        * df['Trib_Rn'] # Radon concentration of inflowing tributary water if present
        * df['dx'] # Reach length
        '''
        try:
            self.sfr
        except:
            print("SFR package not activated in current model")
            return
        # End try

        FlRnEC = Radon_EC_simple(df, Ini_cond, Rn_decay=Rn_decay)
        Flow, Rn, EC = FlRnEC.Fl_Rn_EC_simul()
        return Flow, Rn, EC

    def getHeads(self):
        headobj = self.importHeads()
        times = headobj.get_times()
        head = headobj.get_data(totim=times[-1])
        return head

    def getFinalHeads(self, filename):
        headobj = bf.HeadFile(filename)
        times = headobj.get_times()
        head = headobj.get_data(totim=times[-1])
        return head

    def getRivFlux(self, name):
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
            print('Not a river boundary')

        riv_exchange = {}
        for per in range(str_pers):
            riv_exchange[per] = []
            for cell in river:
                (l, r, c) = (cell[0], cell[1], cell[2])
                if str_pers > 1:
                    riv_exchange[per] += [[riv_flux[per][l][r][c], (l, r, c)]]
                else:
                    riv_exchange[per] += [[riv_flux[0][l][r][c], (l, r, c)]]
                # End if
            # End for
        # End for

        return riv_exchange
    # End getRivFlux()

    def getRivFluxNodes(self, nodes):
        try:
            cbbobj = self.importCbb()
        except IndexError as e:
            raise IndexError("""Error occurred reading cell-by-cell file - check if model converged
            Error: {}""".format(e))
        # End try
        riv_flux = cbbobj.get_data(text='RIVER LEAKAGE', full3D=True)
        times = cbbobj.get_times()
        if type(times) == str:
            str_pers = 1
        else:
            str_pers = len(times)
        # End if

        river = nodes
        riv_exchange = {}
        for per in xrange(str_pers):
            riv_exchange[per] = []
            for cell in river:
                (l, r, c) = (cell[0], cell[1], cell[2])
                if str_pers > 1:
                    riv_exchange[per] += [[riv_flux[per][l][r][c], (l, r, c)]]
                else:
                    riv_exchange[per] += [[riv_flux[0][l][r][c], (l, r, c)]]
                # End if
            # End for
        # End for

        riv_exchange = np.array([x[0] for x in riv_exchange[0] if type(x[0]) == np.float32]).sum()

        return riv_exchange
    # End getRivFluxNodes()

    def getAverageDepthToGW(self, mask=None):

        headobj = self.importHeads()
        times = headobj.get_times()
        head = headobj.get_data(totim=times[-1])

        if mask is not None:
            return np.mean(self.top[mask] - head[0][mask])
        else:
            return np.mean(self.top - head[0])
        # End if
    # End getAverageDepthToGW()

    def loop_over_zone(self, array):
        mesh_1 = self.model_data.model_mesh3D[1]
        arr_zoned = [np.full(mesh_1.shape[1:3], np.nan)] * int(np.max(mesh_1))

        for zone in range(int(np.max(mesh_1))):
            temp = np.array([np.full(mesh_1.shape[1:3], np.nan)] * mesh_1.shape[0])

            zone_1 = float(zone + 1)
            for layer in range(mesh_1.shape[0]):
                mesh_1_layer = mesh_1[layer]
                temp[layer][mesh_1_layer == zone_1] = array[layer][mesh_1_layer == zone_1]
            # End for

            masked_temp = np.ma.masked_array(temp, np.isnan(temp))
            arr_zoned[zone] = np.mean(masked_temp, axis=0)
        # End for

        return arr_zoned
    # End loop_over_zone()

    def ConcsByZone(self, concs):
        return self.loop_over_zone(concs)
    # End ConcsByZone()

    def HeadsByZone(self, heads):
        return self.loop_over_zone(heads)
    # End HeadsByZone()

    def writeObservations(self):

        # Set model output arrays to None to initialise
        head = None
        sfr_df = None
        stream_options = ['stage', 'depth', 'discharge']
        # Write observation to file
        for obs_set in self.model_data.observations.obs_group.keys():
            obs_type = self.model_data.observations.obs_group[obs_set]['obs_type']
            # Import the required model outputs for processing
            if obs_type == 'head':
                # Check if model outputs have already been imported and if not import
                if not head:
                    headobj = self.importHeads()
                    head = headobj.get_alldata()
            elif obs_type in stream_options:
                try:
                    self.sfr_df
                except:
                    sfr_df = self.importSfrOut()
                # End except
            else:
                continue
            # End if

            obs_df = self.model_data.observations.obs_group[obs_set]['time_series']
            obs_df = obs_df[obs_df['active'] == True]
            sim_map_dict = self.model_data.observations.obs_group[obs_set]['mapped_observations']

            if obs_type in stream_options:
                sfr_location = self.model_data.observations.obs_group[obs_set]['locations']['seg_loc']
                for zone in obs_df['zone'].unique():
                    if len(obs_df['zone'].unique()) == 1:
                        zone_txt = obs_set
                    else:
                        zone_txt = obs_set + zone
                    # End if
                    with open(os.path.join(self.data_folder, 'observations_' + zone_txt + '.txt'), 'w') as f:
                        obs_df_zone = obs_df[obs_df['zone'] == zone]
                        for observation in obs_df_zone.index:
                            interval = int(obs_df_zone['interval'].loc[observation])
                            name = obs_df_zone['name'].loc[observation]
                            seg = sfr_location.loc[name]
                            sfr = sfr_df
                            col_of_interest = obs_type
                            if obs_type == 'discharge':
                                col_of_interest = 'Qout'
                            sim_obs = sfr[(sfr['segment'] == seg) &
                                          (sfr['time'] == interval)][col_of_interest]
                            f.write('%f\n' % sim_obs)
                        # End for
                    # End with
                # End for
            # End if

            if obs_type == 'head':
                for zone in obs_df['zone'].unique():
                    if len(obs_df['zone'].unique()) == 1:
                        zone_txt = 'head'
                    else:
                        zone_txt = zone
                    # End if
                    with open(os.path.join(self.data_folder, 'observations_' + zone_txt + '.txt'), 'w') as f:
                        obs_df_zone = obs_df[obs_df['zone'] == zone]
                        for observation in obs_df_zone.index:
                            interval = int(obs_df_zone['interval'].loc[observation])
                            name = obs_df_zone['name'].loc[observation]
                            (x_cell, y_cell) = self.model_data.mesh2centroid2Dindex[
                                (sim_map_dict[name][1], sim_map_dict[name][2])]
                            (lay, row, col) = [sim_map_dict[name][0],
                                               sim_map_dict[name][1], sim_map_dict[name][2]]

                            sim_heads = [head[interval][lay][row][col]]

                            sim_head = np.mean(sim_heads)
                            f.write('%f\n' % sim_head)
                        # End for
                    # End with
                # End for
            # End if
        # End for
    # End writeObservations()

    def getObservation(self, obs, interval, obs_set):
        # Set model output arrays to None to initialise
        head = None

        obs_type = self.model_data.observations.obs_group[obs_set]['obs_type']
        # Import the required model outputs for processing
        if obs_type == 'head':
            # Check if model outputs have already been imported and if not import
            if not head:
                headobj = self.importHeads()
                head = headobj.get_alldata()
        elif obs_type == 'stage':
            pass
        # End if

        sim_map_dict = self.model_data.observations.obs_group[obs_set]['mapped_observations']
        sim_head = head[interval][sim_map_dict[obs][0]][sim_map_dict[obs][1]][sim_map_dict[obs][2]]
        dtw = self.model_data.model_mesh3D[0][0][
            sim_map_dict[obs][1]][sim_map_dict[obs][2]] - sim_head
        return sim_head, dtw
    # End getObservation()

    def compare_observed_head(self, obs_set, simulated, nper=0):
        mesh_1 = self.model_data.model_mesh3D[1]
        obs_group = self.model_data.observations.obs_group[obs_set]
        locations = obs_group['locations']

        self.comparison = []
        self.comp_zone = []
        self.obs_loc_val_zone = []
        obs_df = obs_group['time_series']
        obs_df = obs_df[obs_df['active'] == True]
        obs_df = obs_df[obs_df['interval'] == nper]
        sim_map_dict = obs_group['mapped_observations']
        # self.model_data.observations.obs_group[obs_set]['time_series']['name']:
        for observation in obs_df['name']:
            idx = obs_df[obs_df['name'] == observation].index.tolist()[0]
            obs = obs_df.get_value(idx, 'value')

            mapped_obs = sim_map_dict[observation]
            sim = simulated[mapped_obs[0]][mapped_obs[1]][mapped_obs[2]]

            self.comparison += [[obs, sim]]
            self.comp_zone += mesh_1[sim] / np.max(mesh_1)

            x_y = (locations['Easting'].loc[observation], locations['Northing'].loc[observation])
            zone = mesh_1[sim]
            self.obs_loc_val_zone += [[obs, x_y, zone, sim]]
        # End for
    # End compare_observed_head()

    def CompareObservedHead(self, obs_set, simulated, nper=0):
        """
        :param obs_set: dict, observation set
        :param simulated: dict, simulated result set
        :param nper: int, number of periods
        """
        warnings.warn("Use of deprecated method `CompareObservedHead`, use `compare_observed_head` instead",
                      DeprecationWarning)
        return self.compare_observed_head(obs_set, simulated, nper)
    # End CompareObservedHead()

    def CompareObserved(self, obs_set, simulated, nper=0):
        self.obs_sim_zone = []
        obs_df = self.model_data.observations.obs_group[obs_set]['time_series']
        obs_df = obs_df[obs_df['active'] == True]
        obs_df = obs_df[obs_df['interval'] == nper]
        sim_map_dict = self.model_data.observations.obs_group[obs_set]['mapped_observations']
        # self.model_data.observations.obs_group[obs_set]['time_series']['name']:
        for observation in obs_df['name']:
            idx = obs_df[obs_df['name'] == observation].index.tolist()[0]
            obs = obs_df.get_value(idx, 'value')
            sim = simulated[sim_map_dict[observation][0]][
                sim_map_dict[observation][1]][sim_map_dict[observation][2]]
            zone = self.model_data.model_mesh3D[1][sim_map_dict[observation][0]][
                sim_map_dict[observation][1]][sim_map_dict[observation][2]]
            x = self.model_data.observations.obs_group[obs_set][
                'locations']['Easting'].loc[observation]
            y = self.model_data.observations.obs_group[obs_set][
                'locations']['Northing'].loc[observation]
            if np.isnan(sim):
                print sim, obs, zone
                continue
            self.obs_sim_zone += [[obs, sim, zone, x, y]]
    # End ComparedObserved()

    def importHeads(self, path=None, name=None):
        if path:
            headobj = bf.HeadFile(path + name + '.hds')  # , precision='double')
            return headobj
        else:
            self.headobj = bf.HeadFile(os.path.join(self.data_folder, self.name + '.hds'))
            return self.headobj
        # End if
    # End importHeads()

    def importSfrOut(self, path=None, name=None, ext='.sfr.out'):
        if path:
            sfrout = SfrFile(os.path.join(path, name + ext))
            self.sfr_df = sfrout.get_dataframe()
        else:
            sfrout = SfrFile(os.path.join(self.data_folder, self.name + ext))
            self.sfr_df = sfrout.get_dataframe()
        # End if
        return self.sfr_df
    # End importSfrOut()

    def importCbb(self):
        """Retrieve data in cell-by-cell budget file"""
        self.cbbobj = bf.CellBudgetFile(os.path.join(self.data_folder, self.name + ".cbc"))
        return self.cbbobj
    # End importCbb()

    def waterBalance(self, iter_num, plot=True, save=False, nper=0):

        cbbobj = bf.CellBudgetFile(os.path.join(self.data_folder, self.name + '.cbc'))

        water_balance_components = cbbobj.textlist
        water_balance = {}
        water_balance_summary = {}
        water_bal_summed_titles = []
        water_bal_summed_values = []
        for component in water_balance_components:
            if 'FLOW' in component:
                continue

            component_stripped = component.strip()
            water_balance[component_stripped] = cbbobj.get_data(text=component, full3D=True)  # [-1]

            wb_element = water_balance[component_stripped][nper]
            pos = component_stripped + '_pos'
            neg = component_stripped + '_neg'

            water_balance_summary[pos] = np.sum(wb_element[wb_element > 0])
            water_balance_summary[neg] = np.sum(wb_element[wb_element < 0])

            if water_balance_summary[pos] is np.ma.masked:  # core.MaskedConstant():
                water_balance_summary[pos] = 0.0
            if water_balance_summary[neg] is np.ma.masked:
                water_balance_summary[neg] = 0.0

            water_bal_summed_titles += pos, neg
            water_bal_summed_values += water_balance_summary[pos], water_balance_summary[neg]

        water_bal_summed_titles += ['Error']
        Error = sum(water_bal_summed_values)
        water_bal_summed_values += [Error]

        wat_bal_df = pd.DataFrame(water_bal_summed_values, water_bal_summed_titles)
        wat_bal_df.columns = ['Flux m^3/d']
        wat_bal_df = wat_bal_df[wat_bal_df['Flux m^3/d'] != 0.]

        if plot:
            self.water_balance_plot(iter_num, wat_bal_df, save)
        else:
            return wat_bal_df
        # end if
    # End waterBalance()

    def waterBalanceTS(self, plot=True):

        cbbobj = bf.CellBudgetFile(os.path.join(self.data_folder, self.name + '.cbc'))

        water_balance_components = cbbobj.textlist
        times = cbbobj.get_times()
        water_balance = {}
        water_balance_summary = {}
        water_bal_summed_titles = []
        water_bal_summed_values = []
        for component in water_balance_components:
            if 'FLOW' in component:
                continue

            component_stripped = component.strip()
            water_balance[component_stripped] = cbbobj.get_data(text=component, full3D=True)
            pos = component_stripped + '_pos'
            neg = component_stripped + '_neg'

            water_balance_summary[pos] = []
            water_balance_summary[neg] = []

            for nper in range(len(times)):
                wb_element = water_balance[component_stripped][nper]

                if not np.sum(wb_element[wb_element > 0]) is np.ma.masked:
                    water_balance_summary[pos] += [np.sum(wb_element[wb_element > 0])]
                else:
                    water_balance_summary[pos] += [0.]

                if not np.sum(wb_element[wb_element < 0]) is np.ma.masked:
                    water_balance_summary[neg] += [np.sum(wb_element[wb_element < 0])]
                else:
                    water_balance_summary[neg] += [0.]

            water_bal_summed_titles += pos, neg
            water_bal_summed_values += np.sum(water_balance_summary[pos]), np.sum(water_balance_summary[neg])

        water_bal_summed_titles += ['Error']
        Error = sum(water_bal_summed_values)
        water_bal_summed_values += [Error]

        wat_bal_ts_df = pd.DataFrame(water_balance_summary)
        wat_bal_ts_df = wat_bal_ts_df[[x for x in wat_bal_ts_df.sum().index if wat_bal_ts_df.sum()[x] != 0]]
        wat_bal_ts_df['time'] = pd.to_datetime(self.start_datetime)
        wat_bal_ts_df['time'] = wat_bal_ts_df['time'] + pd.to_timedelta(times, unit='d')

        if plot:
            ax = wat_bal_ts_df.plot(x='time')
            return wat_bal_ts_df, ax
        else:
            return wat_bal_ts_df
        # End if
    # End waterBalanceTS()

# End ModflowModel()


if __name__ == '__main__':
    pass

# End main
