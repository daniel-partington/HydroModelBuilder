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
# allow import from this module to maintain backwards compatibility
from MT3DModel import MT3DModel
from MT3DPostProcess import MT3DPostProcess
from Radon_EC_simple import Radon_EC_simple


class ModflowModel(object):

    def __init__(self, model_data, data_folder=None, **kwargs):
        """

        :param model_data: ModelManager instance, containing all the data for the model
        :param data_folder: str, path to data folder (MODFLOW model run files)
        """
        self.model_data = model_data
        self.name = model_data.name
        if not data_folder:
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
        """Create and setup MODFLOW DIS package.

        See [Flopy MF DIS]_ documentation

        .. [Flopy MF DIS] MODFLOW DIS package
           (https://modflowpy.github.io/flopydoc/mfdis.html)
        """
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
        """Create and setup MODFLOW BAS package.

        See [Flopy MF BAS]_ documentation

        .. [Flopy MF BAS] MODFLOW BAS package
           (https://modflowpy.github.io/flopydoc/mfbas.html)

        :param nlay: int, number of layers
        :param nrow: int, number of rows
        :param ncol: int, number of columns
        """
        ibound = np.copy(self.model_data.model_mesh3D[1])  # Variables for the BAS package
        ibound[ibound == -1] = 0

        self.bas = flopy.modflow.ModflowBas(self.mf, ibound=ibound, strt=self.strt)
        if self.verbose and self.check:
            self.bas.check()
        # End if
    # End setupBASPackage()

    def setupNWTpackage(self, headtol, fluxtol):
        """Create and setup the NWT package.
        See [Flopy MF NWT]_ documentation

        .. [Flopy MF NWT] MODFLOW NWT package
           (https://modflowpy.github.io/flopydoc/mfnwt.html)

        :param headtol: float, maximum head change
        :param fluxtol: float, maximum flux
        """
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
        """Set up the upstream weighting package for the MODFLOW NWT solver.
        See [Flopy MF PCG]_ documentation

        .. [Flopy MF PCG] MODFLOW River package
           (https://modflowpy.github.io/flopydoc/mfpcg.html)
        """
        self.pcg = flopy.modflow.ModflowPcg(self.mf)  # if using mf 2005
    # End setupPCGpachage()

    def setupUPWpackage(self):
        """Set up the upstream weighting package for the MODFLOW NWT solver.
        See [Flopy MF UPW]_ documentation

        .. [Flopy MF UPW] MODFLOW River package
           (https://modflowpy.github.io/flopydoc/mfupw.html)
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
        """Create the MODFLOW river (RIV) package.
        See [Flopy MF Riv]_ documentation

        .. [Flopy MF Riv] MODFLOW River package
           (https://modflowpy.github.io/flopydoc/mfriv.html)

        :param lrcd: dict or None, stress period data. (Default value = None)
        """
        self.riv = flopy.modflow.ModflowRiv(self.mf, ipakcb=53, stress_period_data=lrcd)
        if self.verbose and self.check:
            self.riv.check()
        # End if
    # End createRIVpackage()

    def createSTRpackage(self, STRvariables=None):
        """Create the MODFLOW stream (STR) package.
        See [Flopy MF STR]_ documentation

        .. [Flopy MF STR] MODFLOW stream package
           (https://modflowpy.github.io/flopydoc/mfstr.html)

        :param STRvariables: dict or None, stress period data. (Default value = None)
        """
        self.str = flopy.modflow.ModflowStr(self.mf, ipakcb=53, stress_period_data=STRvariables)
        if self.verbose and self.check:
            self.str.check()
    # End createSTRpackage()

    def createSFRpackage(self, reach_data, segment_data):
        """Create the MODFLOW (SFR2) package.
        See [Flopy MF SFR2]_ documentation

        .. [Flopy MF SFR2] MODFLOW stream-flow routing package
           (https://modflowpy.github.io/flopydoc/mfsfr2.html)

        :param reach_data: ndarray, array holding reach data.
        :param segment_data: ndarray, array holding stream segment data.
        """
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
        """Creates and adds the MODFLOW Gage package.

        See [Flopy mfgage]_ documentation

        .. [Flopy mfgage] MODFLOW Gage
           (https://modflowpy.github.io/flopydoc/mfgage.html)

        :param gages: list or ndarray, data for each gaging location
        :param files:  (Default value = None)
        """
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
        """Creates and adds the MODFLOW Drain package.

        See [Flopy mfdrn]_ documentation

        .. [Flopy mfdrn] MODFLOW DRN
           (https://modflowpy.github.io/flopydoc/mfdrn.html)

        :param lrcsc: list of boundaries, recarrays, or dictionary of boundaries.
                      (Default value = None)
        """
        self.drn = flopy.modflow.ModflowDrn(self.mf, ipakcb=53, stress_period_data=lrcsc)
        if self.verbose and self.check:
            self.drn.check()
    # End createDRNpackage()

    def createGHBpackage(self, lrcsc=None):
        """Creates General Head Boundary module to the model.

        See [Flopy mfghb]_ documentation

        .. [Flopy mfghb] MODFLOW GHB
           (https://modflowpy.github.io/flopydoc/mfghb.html)

        :param lrcsc: list of boundaries, recarray of boundaries or, dictionary of boundaries.
                      (Default value = None)
        """
        self.ghb = flopy.modflow.ModflowGhb(self.mf, ipakcb=53, stress_period_data=lrcsc)
        if self.verbose and self.check:
            self.ghb.check()
    # end createGHBpackage()

    def createRCHpackage(self, rchrate=None):
        """Add RCH package to the MODFLOW model to represent recharge.

        See [Flopy mfrch]_ documentation

        .. [Flopy mfrch] MODFLOW Recharge
           (https://modflowpy.github.io/flopydoc/mfrch.html)

        :param rchrate: float or array of floats, recharge flux (Default value = None)
        """
        self.rch = flopy.modflow.ModflowRch(self.mf, ipakcb=53, rech=rchrate, nrchop=3)
        if self.verbose and self.check:
            self.rch.check()
    # End createRCHpackage()

    def createWELpackage(self, lrcq=None):
        """Add WEL package to the MODFLOW model to represent pumping wells
        Expects a dictionary with an array of well location lay row col and flux
        at each stress period, e.g.

        >>> lrcq = {}
        >>> lrcq[0] = [[0, 7, 7, -100.]] # layer, row, column, flux

        See [Flopy mfwel]_ documentation

        .. [Flopy mfwel] MODFLOW WEL package
           (https://modflowpy.github.io/flopydoc/mfwel.html)

        :param lrcq: dict, stress period data. (Default value = None)
        """
        self.wel = flopy.modflow.ModflowWel(self.mf, ipakcb=53, stress_period_data=lrcq)
        if self.verbose and self.check:
            self.wel.check()
    # End createWElpackage()

    def createOCpackage(self):
        """Add OC (output control) package to the MODFLOW model

        See [Flopy mfoc]_ documentation

        .. [Flopy mfoc] MODFLOW Output Control Module
           (https://modflowpy.github.io/flopydoc/mfoc.html)
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
        """Add LMT package to the MODFLOW model to allow linking with MT3DMS

        See [Flopy mflmt]_ documentation

        .. [Flopy mflmt] MODFLOW Link-MT3DMS
           (https://modflowpy.github.io/flopydoc/mflmt.html)
        """
        self.lmt = flopy.modflow.ModflowLmt(self.mf,
                                            output_file_header='extended',
                                            output_file_format='formatted',
                                            package_flows=['sfr'])
    # End createLMTpackage()

    def finaliseModel(self):
        """Deprecated method"""
        # Write the MODFLOW model input files
        warnings.warn("Deprecated method called. Use `finalize_model()` instead", DeprecationWarning)
        self.finalize_model()
    # end finaliseModel

    def finalize_model(self):
        """Write out inputs used."""
        self.mf.write_input()
    # End finalize_model()

    def add_bc_to_target(self, bc_array, target):
        """TODO: Docs:

        :param bc_array:
        :param target: dict, dictionary to update

        :returns: dict, updated dictionary
        """
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
        """Build MODFLOW model.

        * Creates model discretization (MODFLOW dis package)
        * Sets up the BAS, NWT, UPW, and OC packages
        * Additionally sets up the RCH, DRN, GHB, and SFR packages based on boundary condition properties.

        If necessary, creates linkages to the RIV, WEL, and LMT packages.

        See the relevant docstring for the following methods:
        * :func:`createDiscretisation <flopyInterface.ModflowModel.createDiscretisation>`
        * :func:`createBASpackage <flopyInterface.ModflowModel.createBASpackage>`
        * :func:`createNWTpackage <flopyInterface.ModflowModel.createNWTpackage>`
        * :func:`createUPWpackage <flopyInterface.ModflowModel.createUPWpackage>`
        * :func:`createOCpackage <flopyInterface.ModflowModel.createOCpackage>`
        * :func:`createRCHpackage <flopyInterface.ModflowModel.createRCHpackage>`
        * :func:`createDRNpackage <flopyInterface.ModflowModel.createDRNpackage>`
        * :func:`createGHBpackage <flopyInterface.ModflowModel.createGHBpackage>`
        * :func:`createSFRpackage <flopyInterface.ModflowModel.createSFRpackage>`
        * :func:`createRIVpackage <flopyInterface.ModflowModel.createRIVpackage>`
        * :func:`createWELpackage <flopyInterface.ModflowModel.createWELpackage>`
        * :func:`createLMTpackage <flopyInterface.ModflowModel.createLMTpackage>`

        :param transport: bool, use/setup the MODFLOW transport package  (Default value = False)
        :param write: bool, writes out the inputs generated (Default value = True)
        :param verbose: bool, print out verbose messages (Default value = True)
        :param check: bool, check MODFLOW output for errors/convergence (Default value = False)
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

        river = {}
        wel = {}

        bc = self.model_data.boundaries.bc
        for boundary in bc:
            bc_boundary = bc[boundary]
            bc_type = bc_boundary['bc_type']
            bc_array = bc_boundary['bc_array']

            if bc_type == 'recharge':
                self.createRCHpackage(rchrate=bc_array)
            elif bc_type == 'drain':
                self.createDRNpackage(bc_array)
            elif bc_type == 'general head':
                self.createGHBpackage(bc_array)
            elif (bc_type == 'river') or (bc_type == 'channel'):
                river = self.add_bc_to_target(bc_array, river)
            elif (bc_type == 'river_flow'):
                self.createSFRpackage(bc_array[0], bc_array[1])
            elif bc_type == 'wells':
                wel = self.add_bc_to_target(bc_array, wel)
            # End if

        # End for

        if river:
            self.createRIVpackage(river)

        if wel:
            self.createWELpackage(wel)

        if transport:
            self.createLMTpackage()

        if write:
            self.finalize_model()

        if self.verbose and self.check:
            self.checkMODFLOW()
        # End if

    # End buildMODFLOW()

    def checkMODFLOW(self):
        """Check model data for common errors."""
        self.mf.check()
    # End checkMODFLOW

    def runMODFLOW(self, silent=True):
        """Run the MODFLOW model.
        Calls the `checkConvergence` method to capture model run failures that are
        not picked up by flopy.

        :param silent: bool, optional argument for suppressing any output to screen. (Default value = True)

        :returns: bool, successful run with convergence.
        """
        success, buff = self.mf.run_model(silent=silent)
        return self.checkConvergence(fail=not success)
    # End runMODFLOW()

    def checkConvergence(self, path=None, name=None, fail=False):
        """TODO: Docs

        :param path:  (Default value = None)
        :param name:  (Default value = None)
        :param fail:  (Default value = False)

        :returns: bool, model converged (`True`) or failed (`False`)
        """
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

    def Calculate_Rn_from_SFR_with_simple_model(self, df, ini_cond, Rn_decay=0.181):
        """Use a simple model to calculate Radon concentrations in the stream
        based on outputs from the SFR package and using some other data that
        is required as arguments to this function.

        In particular the df is a pandas dataframe that should contain the output
        from the sfr package which can be imported via the sfroutputfile util
        from flopy. The other parameters required are:

        * Ini_Cond = ini_cond  # 3 item list containing Initial flow, radon and ec concentrations
        * Rn_decay = Rn_decay  # constant for Radon decay

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

        :param df: Pandas DataFrame, output from the SFR package
        :param Ini_cond: list, [intial flow, radon, and EC concentration]
        :param Rn_decay: float, constant for radon decay. (Default value = 0.181)

        :returns: tuple, (flow, radon, EC)
        """

        try:
            self.sfr
        except Exception:
            print("SFR package not activated in current model")
            return
        # End try

        FlRnEC = Radon_EC_simple(df, ini_cond, Rn_decay=Rn_decay)
        Flow, Rn, EC = FlRnEC.Fl_Rn_EC_simul()

        return Flow, Rn, EC
    # End Calculate_Rn_from_SFR_with_simple_model()

    def get_heads(self, headobj=None):
        """Get latest head level.

        :param headobj: flopy head object, (Default value = None)

        :returns: ndarray, head levels
        """
        if not headobj:
            headobj = self.import_heads()
        times = headobj.get_times()
        return headobj.get_data(totim=times[-1])
    # End get_heads()

    def get_final_heads(self, filename):
        """Get final head (groundwater levels) at end of simulation from given file.

        :param filename: str, name of head file.

        :returns: ndarray, head levels
        """
        headobj = bf.HeadFile(filename)
        return self.get_heads(headobj)
    # End get_final_heads()

    def getHeads(self, headobj=None):
        """Deprecated method.

        :param headobj: Default value = None)

        :returns: ndarray, head level
        """
        warnings.warn("Use of deprecated method `getHeads`, use `get_heads` instead",
                      DeprecationWarning)
        return self.get_heads(headobj)
    # End getHeads()

    def getFinalHeads(self, filename):
        """Deprecated method.
        :param filename: returns: ndarray, head level
        """
        warnings.warn("Use of deprecated method `getFinalHeads`, use `get_final_heads` instead",
                      DeprecationWarning)
        return self.get_final_heads(filename)
    # End getFinalHeads()

    def getRivFlux(self, name):
        """Get river flux data.

        :param name: str, name of river boundary

        :returns: dict, river exchange values
        """
        cbbobj = self.importCbb()
        riv_flux = cbbobj.get_data(text='RIVER LEAKAGE', full3D=True)

        times = cbbobj.get_times()
        if type(times) == str:
            str_pers = 1
        else:
            str_pers = len(times)
        # End if

        if self.model_data.boundaries.bc[name]['bc_type'] == 'river':
            time_key = self.model_data.boundaries.bc[name]['bc_array'].keys()[0]
            river = self.model_data.boundaries.bc[name]['bc_array'][time_key]
        else:
            print('Not a river boundary')
        # End if

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
        """Get river flux nodes.

        :param nodes:

        :returns: dict, river exchange nodes
        """
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

        riv_exchange = np.array([x[0] for x in riv_exchange[0] if type(x[0]) == np.float32]).sum()

        return riv_exchange
    # End getRivFluxNodes()

    def get_average_depth_to_GW(self, mask=None):
        """Get average head values.

        :param mask:  (Default value = None)

        :returns: float, average head value
        """
        head = self.get_heads()
        if mask:
            return np.mean(self.top[mask] - head[0][mask])
        else:
            return np.mean(self.top - head[0])
        # End if
    # End get_average_depth_to_GW()

    def getAverageDepthToGW(self, mask=None):
        """Deprecated method.

        :param mask: (Default value = None)
        """
        warnings.warn("Use of deprecated method `getAverageDepthToGW`, use `get_average_depth_to_GW` instead",
                      DeprecationWarning)
        return self.get_average_depth_to_GW(mask)
    # End getAverageDepthToGW()

    def loop_over_zone(self, array):
        """Generate a masked array and retrieve average values.

        :param array: ndarray, representing zone.

        :returns: ndarray, average values for zone.
        """
        mesh_1 = self.model_data.model_mesh3D[1]
        rows = mesh_1.shape[0]
        arr_zoned = [np.full(mesh_1.shape[1:3], np.nan)] * int(np.max(mesh_1))

        for zone in range(int(np.max(mesh_1))):
            temp = np.array([np.full(mesh_1.shape[1:3], np.nan)] * rows)

            zone_1 = float(zone + 1)
            for layer in range(rows):
                mesh_1_layer = mesh_1[layer]
                temp[layer][mesh_1_layer == zone_1] = array[layer][mesh_1_layer == zone_1]
            # End for

            masked_temp = np.ma.masked_array(temp, np.isnan(temp))
            arr_zoned[zone] = np.mean(masked_temp, axis=0)
        # End for

        return arr_zoned
    # End loop_over_zone()

    def concs_by_zone(self, concs):
        """TODO: Docs

        :param concs:
        """
        return self.loop_over_zone(concs)
    # End concs_by_zone()

    def ConcsByZone(self, concs):
        """Deprecated method.
        :param concs:

        :returns: ndarray
        """
        warnings.warn("Use of deprecated method `ConcsByZone`, use `concs_by_zone` instead",
                      DeprecationWarning)
        return self.concs_by_zone(concs)
    # End ConcsByZone()

    def heads_by_zone(self, heads):
        """Retrieve average head values for each zone.

        :param heads:

        :returns: ndarray
        """
        return self.loop_over_zone(heads)
    # End heads_by_zone()

    def HeadsByZone(self, heads):
        """Deprecated method.
        :param heads:
        """
        warnings.warn("Use of deprecated method `HeadsByZone`, use `heads_by_zone` instead",
                      DeprecationWarning)
        return self.heads_by_zone(heads)
    # End HeadsByZone()

    def writeObservations(self):
        """Write out observation data."""
        # Set model output arrays to None to initialise
        head = None
        sfr_df = None
        stream_options = ['stage', 'depth', 'discharge']
        # Write observation to file
        obs_group = self.model_data.observations.obs_group
        for obs_set in obs_group:
            obs_type = obs_group[obs_set]['obs_type']
            # Import the required model outputs for processing
            if obs_type == 'head':
                # Check if model outputs have already been imported and if not import
                if not head:
                    headobj = self.import_heads()
                    head = headobj.get_alldata()
            elif obs_type in stream_options:
                try:
                    self.sfr_df
                except Exception:
                    sfr_df = self.importSfrOut()
                # End except
            else:
                continue
            # End if

            group_set = obs_group[obs_set]
            obs_df = self.get_active_obs_group(group_set)
            sim_map_dict = group_set['mapped_observations']

            unique_obs_df_zone = obs_df['zone'].unique()
            if obs_type in stream_options:
                sfr_location = group_set['locations']['seg_loc']
                for zone in unique_obs_df_zone:
                    zone_txt = obs_set if len(unique_obs_df_zone) == 1 else obs_set + zone
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

                            sim_obs = sfr[(sfr['segment'] == seg) & (sfr['time'] == interval)][col_of_interest]
                            f.write('%f\n' % sim_obs)
                        # End for
                    # End with
                # End for
            elif obs_type == 'head':
                for zone in unique_obs_df_zone:
                    zone_txt = 'head' if len(unique_obs_df_zone) == 1 else zone
                    with open(os.path.join(self.data_folder, 'observations_' + zone_txt + '.txt'), 'w') as f:
                        obs_df_zone = obs_df[obs_df['zone'] == zone]
                        for observation in obs_df_zone.index:
                            interval = int(obs_df_zone['interval'].loc[observation])
                            name = obs_df_zone['name'].loc[observation]
                            name_dict = sim_map_dict[name]
                            (x_cell, y_cell) = self.model_data.mesh2centroid2Dindex[(name_dict[1], name_dict[2])]
                            (lay, row, col) = [name_dict[0], name_dict[1], name_dict[2]]

                            sim_heads = [head[interval][lay][row][col]]
                            sim_head = np.mean(sim_heads)
                            f.write('%f\n' % sim_head)
                        # End for
                    # End with
                # End for
            else:
                print("Unknown observation type!")
            # End if
        # End for
    # End writeObservations()

    def get_observation(self, obs, interval, obs_set):
        """Get observation data.

        :param obs:
        :param interval:
        :param obs_set:

        :returns: tuple,
        """
        # Set model output arrays to None to initialise
        head = None

        obs_group = self.model_data.observations.obs_group[obs_set]

        obs_type = obs_group['obs_type']
        # Import the required model outputs for processing
        if obs_type == 'head':
            # Check if model outputs have already been imported and if not import
            if not head:
                headobj = self.import_heads()
                head = headobj.get_alldata()
        elif obs_type == 'stage':
            pass
        # End if

        sim_map_dict = obs_group['mapped_observations']
        sim_obs = sim_map_dict[obs]
        sim_head = head[interval][sim_obs[0]][sim_obs[1]][sim_obs[2]]
        dtw = self.model_data.model_mesh3D[0][0][sim_obs[1]][sim_obs[2]] - sim_head
        return sim_head, dtw
    # End get_observation()

    def getObservation(self, obs, interval, obs_set):
        """Deprecated method.

        :param obs:
        :param interval:
        :param obs_set:
        """
        warnings.warn("Use of deprecated method `getObservation`, use `get_observation` instead",
                      DeprecationWarning)
        return self.get_observation(obs, interval, obs_set)
    # End getObservation()

    def get_active_obs_group(self, obs_group, nper=None):
        """Get active observations within a group.

        :param obs_group: dict, holding observation time series
        :param nper: int, optional period number (Default value = None)

        :returns: DataFrame
        """
        assert ('time_series' in obs_group) and \
               ('active' in obs_group['time_series']) and \
               ('interval' in obs_group['time_series']), \
            "Given DataFrame is missing required columns"

        obs_df = obs_group['time_series']
        obs_df = obs_df[obs_df['active'] == True]

        if nper:
            obs_df = obs_df[obs_df['interval'] == nper]
        # End if

        return obs_df
    # End get_active_obs_group()

    def compare_observed_head(self, obs_set, simulated, nper=0):
        """Create and set data for groundwater head comparison.

        :param obs_set: str, name of observation set
        :param simulated: ndarray, simulated observation set
        :param nper: int, period number. (Default value = 0)

        :returns: list[list], data as assigned to `obs_loc_val_zone`
        """
        mesh_1 = self.model_data.model_mesh3D[1]
        obs_group = self.model_data.observations.obs_group[obs_set]
        locations = obs_group['locations']

        self.comparison = []
        self.comp_zone = []
        self.obs_loc_val_zone = []

        obs_df = self.get_active_obs_group(obs_group, nper)
        sim_map_dict = obs_group['mapped_observations']
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

        return self.obs_loc_val_zone
    # End compare_observed_head()

    def CompareObservedHead(self, obs_set, simulated, nper=0):
        """ Deprecated method.

        :param obs_set:
        :param simulated:
        :param nper:  (Default value = 0)
        """
        warnings.warn("Use of deprecated method `CompareObservedHead`, use `compare_observed_head` instead",
                      DeprecationWarning)
        return self.compare_observed_head(obs_set, simulated, nper)
    # End CompareObservedHead()

    def compare_observed(self, obs_set, simulated, nper=0):
        """Create a combined observed and simulated dataset and assigns to the `obs_sim_zone` attribute.

        :param obs_set: str, name of observation set
        :param simulated: ndarray, simulated observation set
        :param nper: int, period number. (Default value = 0)

        :returns: list[list], as set in `obs_sim_zone`
        """
        self.obs_sim_zone = []
        obs_group = self.model_data.observations.obs_group[obs_set]
        obs_df = self.get_active_obs_group(obs_group, nper)
        sim_map_dict = obs_group['mapped_observations']

        for observation in obs_df['name']:
            idx = obs_df[obs_df['name'] == observation].index.tolist()[0]
            obs = obs_df.get_value(idx, 'value')

            sim_obs = sim_map_dict[observation]
            sim = simulated[sim_obs[0]][sim_obs[1]][sim_obs[2]]
            zone = self.model_data.model_mesh3D[1][sim_obs[0]][sim_obs[1]][sim_obs[2]]
            x = obs_group['locations']['Easting'].loc[observation]
            y = obs_group['locations']['Northing'].loc[observation]
            if np.isnan(sim):
                print("Sim is NaN:", sim, obs, zone, "|", observation)
                continue
            self.obs_sim_zone += [[obs, sim, zone, x, y]]
        # End for

        return self.obs_sim_zone
    # End compare_observed()

    def CompareObserved(self, obs_set, simulated, nper=0):
        """ Deprecated method.

        :param obs_set:
        :param simulated:
        :param nper:  (Default value = 0)
        """
        warnings.warn("Use of deprecated method `CompareObserved`, use `compare_observed` instead",
                      DeprecationWarning)
        return self.compare_observed(obs_set, simulated, nper)
    # End ComparedObserved()

    def import_heads_from_file(self, path=None, name=None):
        """Import head data from specified file.

        :param path: str, path to file. (Default value = None)
        :param name: str, filename to load data from. (Default value = None)

        :returns: flopy headobj data
        """
        if path:
            headobj = bf.HeadFile(os.path.join(path, name + '.hds'))
        else:
            headobj = self.import_heads()
        # End if

        return headobj
    # End import_heads_from_file()

    def import_heads(self):
        """Import head data associated with this model.

        :returns: flopy headobj data
        """
        if not hasattr(self, 'headobj'):
            self.headobj = bf.HeadFile(os.path.join(self.data_folder, self.name + '.hds'))
        # End if

        return self.headobj
    # End import_heads()

    def importHeads(self, path=None, name=None):
        """ Deprecated method.
        :param path: Default value = None)
        :param name:  (Default value = None)
        :param name:  (Default value = None)
        """
        warnings.warn("""Use of method that will be removed in the future.
                      Use `import_heads()` or `import_heads_from_file()`""", FutureWarning)
        if not path:
            warnings.warn("Deprecated method called. Use `import_heads()` instead", DeprecationWarning)
            return self.import_heads()
        else:
            warnings.warn("Deprecated method called. Use `import_heads_from_file()` instead", DeprecationWarning)
            return self.import_heads_from_file(path, name)
    # End importHeads()

    def importSfrOut(self, path=None, name=None, ext='.sfr.out'):
        """Deprecated method.

        :param path:  (Default value = None)
        :param name:  (Default value = None)
        :param ext:  (Default value = '.sfr.out')
        """
        warnings.warn("Deprecated method called. Use `import_sfr_out()` instead", DeprecationWarning)
        return self.import_sfr_out(path, name, ext)
    # End importSfrOut()

    def import_sfr_out(self, path=None, name=None, ext='.sfr.out'):
        """TODO: Docs

        :param path:  (Default value = None)
        :param name:  (Default value = None)
        :param ext:  (Default value = '.sfr.out')

        :returns: DataFrame, SFR data
        """
        if path:
            sfrout = SfrFile(os.path.join(path, name + ext))
            self.sfr_df = sfrout.get_dataframe()
        else:
            sfrout = SfrFile(os.path.join(self.data_folder, self.name + ext))
            self.sfr_df = sfrout.get_dataframe()
        # End if

        return self.sfr_df
    # End import_sfr_out()

    def import_cbb(self):
        """Retrieve data in cell-by-cell budget file"""

        self.cbbobj = bf.CellBudgetFile(os.path.join(self.data_folder, self.name + ".cbc"))
        return self.cbbobj
    # End import_cbb()

    def importCbb(self):
        """Deprecated method."""
        warnings.warn("Use of deprecated method `importCbb`, use `import_cbb` instead",
                      DeprecationWarning)
        return self.import_cbb()
    # End importCbb()

    def waterBalance(self, iter_num, plot=True, save=False, nper=0):
        """TODO: Docs

        :param iter_num:
        :param plot:  (Default value = True)
        :param save:  (Default value = False)
        :param nper:  (Default value = 0)
        """
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
        wat_bal_df = wat_bal_df[wat_bal_df['Flux m^3/d'] != 0.0]

        if plot:
            self.water_balance_plot(iter_num, wat_bal_df, save)
        else:
            return wat_bal_df
        # End if
    # End waterBalance()

    def waterBalanceTS(self, plot=True):
        """TODO Docs

        :param plot:  (Default value = True)
        """

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
