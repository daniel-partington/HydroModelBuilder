import inspect
import os
import warnings
from types import MethodType

# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# import flopy
import flopy.utils.binaryfile as bf
import viz.MT3D_PP_viz as pp_viz  # Visualization extension methods


class MT3DPostProcess(object):

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
        self.concobj = bf.UcnFile(os.path.join(self.mf_model.data_folder, 'MT3D001.UCN'))
        return self.concobj
    # End importConcs()

    def importSftConcs(self):
        self.sft_conc = pd.read_csv(os.path.join(
            self.mf_model.data_folder,
            self.mt_name),
            delim_whitespace=True,
            skiprows=1)
        return self.sft_conc
    # End importSftConcs()

    def ConcsByZone(self, concs):

        self.mf_model.model_data.model_mesh3D[1]

        concs_zoned = [np.full(self.mf_model.model_data.model_mesh3D[1].shape[
                               1:3], np.nan)] * int(np.max(self.mf_model.model_data.model_mesh3D[1]))

        for zone in range(int(np.max(self.mf_model.model_data.model_mesh3D[1]))):
            temp_concs = np.array([np.full(self.mf_model.model_data.model_mesh3D[1].shape[
                                  1:3], np.nan)] * self.mf_model.model_data.model_mesh3D[1].shape[0])
            # for each layer in mesh get the heads from zone and average:
            for layer in range(self.mf_model.model_data.model_mesh3D[1].shape[0]):
                temp_concs[layer][self.mf_model.model_data.model_mesh3D[1][layer] == float(
                    zone + 1)] = concs[layer][self.mf_model.model_data.model_mesh3D[1][layer] == float(zone + 1)]
            masked_temp_concs = np.ma.masked_array(temp_concs, np.isnan(temp_concs))
            concs_zoned[zone] = np.mean(masked_temp_concs, axis=0)
        # End for
        return concs_zoned
    # End ConcsByZone()

    def compare_observed(self, obs_set, simulated, nper=0):
        self.obs_sim_zone = []
        obs_df = self.mf_model.model_data.observations.obs_group[obs_set]['time_series']
        obs_df = obs_df[obs_df['active'] == True]
        # obs_df = obs_df[obs_df['interval'] == nper]
        sim_map_dict = self.mf_model.model_data.observations.obs_group[
            obs_set]['mapped_observations']
        # self.model_data.observations.obs_group[obs_set]['time_series']['name']:
        for observation in obs_df['name']:
            idx = obs_df[obs_df['name'] == observation].index.tolist()[0]
            obs = obs_df.get_value(idx, 'value')
            sim = simulated[sim_map_dict[observation][0]][
                sim_map_dict[observation][1]][sim_map_dict[observation][2]]
            zone = self.mf_model.model_data.model_mesh3D[1][sim_map_dict[observation][0]][
                sim_map_dict[observation][1]][sim_map_dict[observation][2]]
            x = self.mf_model.model_data.observations.obs_group[
                obs_set]['locations']['Easting'].loc[observation]
            y = self.mf_model.model_data.observations.obs_group[
                obs_set]['locations']['Northing'].loc[observation]
            if np.isnan(sim):
                print(sim, obs, zone)
                continue
            self.obs_sim_zone += [[obs, sim, zone, x, y]]
        # End for
    # End compare_observed()

    def CompareObserved(self, obs_set, simulated, nper=0):
        warnings.warn("Use of deprecated method `CompareObserved`, use `compare_observed` instead",
                      DeprecationWarning)
        return self.compare_observed(obs_set, simulated, nper)
    # End CompareObserved()

    def writeObservations(self, specimen):

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
                        concobj = self.importConcs()
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
