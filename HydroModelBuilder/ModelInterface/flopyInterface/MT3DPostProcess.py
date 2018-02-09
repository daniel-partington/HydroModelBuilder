import inspect
import os
import warnings
from types import MethodType

# import flopy
import flopy.utils.binaryfile as bf
# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import viz.MT3D_PP_viz as pp_viz  # Visualization extension methods
from flopyInterface import ModflowModel


class MT3DPostProcess(ModflowModel):
    """TODO: Docs"""

    def __init__(self, mf_model, mt_name=None):
        self.mf_model = mf_model
        self.mt_name = mt_name

        # Add visualization methods to this instance
        viz_methods = inspect.getmembers(pp_viz, inspect.isfunction)
        for method_name, func in viz_methods:
            setattr(self, method_name, MethodType(func, self))
        # End for

    # End __init__()

    def importConcs(self):
        """TODO: Docs"""
        warnings.warn("Use of deprecated method `importConcs`, use `import_concs` instead",
                      DeprecationWarning)
        return self.import_concs()
    # End importConcs()

    def import_concs(self):
        """TODO: Docs"""
        self.concobj = bf.UcnFile(os.path.join(self.mf_model.data_folder, 'MT3D001.UCN'))
        return self.concobj
    # End import_concs()

    def importSftConcs(self):
        """TODO: Docs"""
        self.sft_conc = pd.read_csv(os.path.join(
            self.mf_model.data_folder,
            self.mt_name),
            delim_whitespace=True,
            skiprows=1)
        return self.sft_conc
    # End importSftConcs()

    def writeObservations(self, specimen):
        """
        :param specimen:
        """

        # Set model output arrays to None to initialise
        conc = None
        sft_conc = None
        data_folder = self.mf_model.data_folder
        model_data = self.mf_model.model_data
        obs_group = self.mf_model.model_data.observations.obs_group

        # Write observation to file
        for obs_set in model_data.observations.obs_group.keys():
            obs_type = obs_group[obs_set]['obs_type']
            # Import the required model outputs for processing
            if obs_type not in ['concentration', 'EC', 'Radon']:
                continue
            else:
                if (obs_type == 'concentration') & (specimen == 'C14'):
                    # Check if model outputs have already been imported and if not import
                    if not conc:
                        concobj = self.import_concs()
                        conc = concobj.get_alldata()  # (totim=times[0])
                    # End if
                elif (obs_type == 'EC') & (specimen == 'EC'):
                    try:
                        sft_conc['TIME']
                    except:
                        sft_conc = self.importSftConcs()
                elif obs_type == 'Radon':
                    continue
                else:
                    continue
                # End if
            # End if

            obs_df = obs_group[obs_set]['time_series']
            obs_df = obs_df[obs_df['active'] == True]
            sim_map_dict = obs_group[obs_set]['mapped_observations']

            with open(data_folder + os.path.sep + 'observations_' + obs_set + '.txt', 'w') as f:
                if obs_group[obs_set]['domain'] == 'stream':
                    sft_location = self.mf_model.model_data.observations.obs_group[obs_set]['locations']['seg_loc']

                    for observation in obs_df.index:
                        interval = int(obs_df['interval'].loc[observation])
                        name = obs_df['name'].loc[observation]
                        seg = sft_location.loc[name]
                        sft = sft_conc
                        times = sft['TIME'].unique()
                        col_of_interest = obs_type
                        if obs_type == 'EC':
                            col_of_interest = 'SFR-CONCENTRATION'
                        sim_conc = sft[(sft['SFR-NODE'] == seg) &
                                       (sft['TIME'] == times[interval])][col_of_interest]
                        f.write('%f\n' % sim_conc)

                if obs_group[obs_set]['domain'] == 'porous':
                    for observation in obs_df.index:
                        interval = int(obs_df['interval'].loc[observation])
                        name = obs_df['name'].loc[observation]

                        (x_cell, y_cell) = model_data.mesh2centroid2Dindex[
                            (sim_map_dict[name][1], sim_map_dict[name][2])]

                        (lay, row, col) = [sim_map_dict[name][0],
                                           sim_map_dict[name][1], sim_map_dict[name][2]]
                        sim_conc = [conc[interval][lay][row][col]]
                        sim_conc = np.mean(sim_conc)
                        f.write('%f\n' % sim_conc)
                    # End for
                # End if
            # End with
        # End for
    # End writeObservations()

# End MT3DPostProcess()
