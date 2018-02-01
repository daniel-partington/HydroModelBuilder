import warnings

import flopy
import numpy as np


class MT3DModel(object):

    def __init__(self, mf_model):
        self.mf_model = mf_model

    def createBTNpackage(self):
        # Add the BTN package to the model
        ibound = self.model_data.model_mesh3D[1]
        ibound[ibound == -1] = 0
        self.btn = flopy.mt3d.Mt3dBtn(self.mt, icbund=ibound, ncomp=1, mcomp=1,
                                      cinact=-9.9E1, thkmin=-1.0E-6, ifmtcn=5,
                                      ifmtnp=0, ifmtrf=0, ifmtdp=0, nprs=0,
                                      timprs=None, savucn=1, nprobs=0,
                                      chkmas=1, nprmas=1, dt0=10000.0, ttsmax=100000.0)

    def createADVpackage(self):
        # Add the ADV package to the model
        self.adv = flopy.mt3d.Mt3dAdv(self.mt, mixelm=0, percel=1,
                                      mxpart=250000, nadvfd=1, itrack=3,
                                      wd=0.5, dceps=1.0E-4, nplane=0, npl=5,
                                      nph=8, npmin=1, npmax=16, nlsink=0,
                                      npsink=8, dchmoc=0.01)

    def createDSPpackage(self):
        # Add the DSP package to the model
        self.dfp = flopy.mt3d.Mt3dDsp(self.mt, multiDiff=True, al=10.,
                                      trpt=0.1, trpv=0.1, dmcoef=0.0)

    def createRCTpackage(self):
        # Add the RCT package to the model
        self.rct = flopy.mt3d.Mt3dRct(self.mt, isothm=1, ireact=1,
                                      igetsc=0, rc1=np.log(2) / (5730 * 365))

    def createGCGpackage(self):
        # Add the GCG package to the model
        self.gcg = flopy.mt3d.Mt3dGcg(self.mt, mxiter=1000, iter1=100, isolve=1,
                                      ncrs=0, accl=1, cclose=1.0E-4, iprgcg=0)

    def finaliseMT3Dmodel(self):
        warnings.warn("Deprecated. Use `finalize_MT3D_model()` instead.", DeprecationWarning)
        self.finalize_MT3D_model()
    # end finaliseMT3Dmodel

    def finalize_MT3D_model(self):
        self.mt.write_input()
    # End finalize_MT3D_model()

    def buildMT3D(self):

        self.mt = flopy.mt3d.Mt3dms(modelname=self.name + '_transport', ftlfilename='mt3d_link.ftl',
                                    modflowmodel=self.mf, model_ws=self.data_folder, exe_name='MT3D-USGS_64.exe')

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
        # End for

        river = {}
        river[0] = []
        bc = self.model_data.boundaries.bc
        for boundary in bc:
            bc_boundary = bc[boundary]
            bc_type = bc_boundary['bc_type']
            bc_array = bc_boundary['bc_array']

            if (bc_type == 'river') or (bc_type == 'channel'):
                time_key = bc_array.keys()[0]
                river[0] += bc_array[time_key]
            # End if

            if bc_type == 'recharge':
                for key in bc_array:
                    self.crch[key] = np.ones_like(bc_array[key])
                    self.crch[key] = self.crch[key] * 100.0
                    self.crch[key][ibound[0] == 0] = 0.0
                # End for
            # End if

            if bc_type == 'wells':
                for key in bc_array:
                    for well in bc_array[key]:
                        ssm_data[key].append((well[0], well[1], well[2], 100.0, itype['WEL']))

            if bc_type == 'general head':
                for key in bc_array:
                    for ghb in bc_array[key]:
                        ssm_data[key].append((ghb[0], ghb[1], ghb[2], 0.0, itype['GHB']))
                    # End for
                # End for
            # End if
        # End for

        if len(river) > 0:
            for key in river:
                for riv in river[key]:
                    ssm_data[key].append((riv[0], riv[1], riv[2], 100.0, itype['RIV']))
                # End for
            # End for
        # End if

        self.ssm = flopy.mt3d.Mt3dSsm(self.mt, stress_period_data=ssm_data, crch=self.crch)
        self.finaliseMT3Dmodel()
    # End buildMT3D()

    def runMT3D(self, silent=False):
        success, buff = self.mt.run_model(silent=silent)
        return success
    # End runMT3D()

# End MT3DModel()
