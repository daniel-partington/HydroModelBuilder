import os

import flopy
import flopy.utils.binaryfile as bf
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import six
from flopy.utils.sfroutputfile import SfrFile
from matplotlib import colors

from HydroModelBuilder.Utilities.model_assessment import (metric_me,
                                                          metric_pbias,
                                                          metric_rmse,
                                                          metric_srms,
                                                          plot_obs_vs_sim)


def cell_bc_to_3D_array_for_plot(self, boundary, bc_array, nper, val_pos):
    '''Convert cell data to 3D array for 3D plotting

    :param boundary: str, name of boundary condition
    :param bc_array: dict of list data for different boundaries including
                     cell location layer, row, column and other pertinent
                     data required by MODFLOW depending on the package
    :param nper: int, time period to use from bc_array
    :param val_pos: int, target value position in bc_array
    '''
    self.bc_3D_array[boundary] = np.zeros_like(self.model_data.model_mesh3D[1])
    bc_cells = [[lrc[0], lrc[1], lrc[2]] for lrc in bc_array[nper]]
    vals = [val[val_pos] for val in bc_array[nper]]
    for index, c in enumerate(bc_cells):
        self.bc_3D_array[boundary][c[0]][c[1]][c[2]] = vals[index]
        self.bc_3D_array['bcs'][c[0]][c[1]][c[2]] = self.bc_counter
    # End for
    self.bc_counter += 1
# End cell_bc_to_3D_array_for_plot()


def create_3D_bc_array(self):
    """Create dict of 3D arrays for visualiation of boundary conditions, which
    can then be passed to the writer.

    """
    mesh3D_1 = self.model_data.model_mesh3D[1]
    mesh_template = np.zeros_like(mesh3D_1)
    self.bc_3D_array = {'zone': np.copy(mesh3D_1)}
    self.bc_3D_array.update({'bcs': mesh_template.copy()})
    self.bc_counter = 1

    bc = self.model_data.boundaries.bc
    for boundary in bc:
        bc_boundary = bc[boundary]
        bc_type = bc_boundary['bc_type']
        bc_array = bc_boundary['bc_array']

        nper = 0
        val_pos = 3

        if bc_type in ['drain', 'general head', 'river', 'channel', 'river_flow', 'wells']:
            if bc_type == 'river_flow':
                val_pos = 5
            elif bc_type == 'wells':
                nper = max(bc_array.keys())
            # End if
            self.cell_bc_to_3D_array_for_plot(boundary, bc_array, nper, val_pos)
        elif bc_type == 'recharge':
            self.bc_3D_array[boundary] = mesh_template.copy()
            if 'zonal_array' in bc_boundary:
                self.bc_3D_array[boundary + '_zonal'] = mesh_template.copy()
                self.bc_3D_array[boundary + '_zonal'][0] = bc_boundary['zonal_array']
                self.bc_3D_array[boundary + '_zonal'][mesh3D_1 == -1] = 0.0
            # End if
            self.bc_3D_array[boundary][0] = bc_array[0]
        # End if
    # End for

# End create_3D_bc_array()


def plotRiverFromRiverSegData(self, ax, names=None, **kwargs):
    """
    :param ax: param names:  (Default value = None)

    :param names:  (Default value = None)

    :param **kwargs:
    """

    river_mapping = self.model_data.river_mapping
    if names is not None:
        keys = [i for i in names if i in river_mapping.keys()]
        not_a_key = [j for j in names if j not in river_mapping.keys()]
        if len(not_a_key) == 0:
            print("Bad names passed that are not in dict: {}".format(not_a_key))
    else:
        keys = river_mapping.keys()
    # End if

    for key in keys:
        river_seg = river_mapping[key]
        all_riv = []
        for x in river_seg['amalg_riv_points_collection']:
            all_riv += x
        # End for
        x = [x[0] for x in all_riv]
        y = [y[1] for y in all_riv]
        ax.plot(x, y, **kwargs)
    # End for
# End plotRiverFromRiverSegData()


def water_balance_plot(self, iter_num, wat_bal_df, save):
    """Plot water balance data.

    :param iter_num: int, iteration number to use in filename

    :param wat_bal_df: DataFrame, water balance data

    :param save: bool, save plot or not
    """

    # Setup params to get water balance aspect ratio looking nice
    # aspect = float(12.5715/((wat_bal_df.max()[0]-wat_bal_df.min()[0])/float(wat_bal_df.shape[1])))

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title('Water Balance')
    wat_bal_df.plot(kind='bar', ax=plt.gca())
    ax.grid(True)
    gridlines = ax.get_xgridlines()
    for line in gridlines:
        line.set_linestyle('-')

    fig.subplots_adjust(left=0.1, right=0.9, bottom=0.35, top=0.95, wspace=0.1, hspace=0.12)
    if save:
        plt.savefig('run_wb_{}.png'.format(iter_num), bbox_inches='tight')
    # End if
# End water_balance_plot()

def SFR_get_data_at_time():
    pass    

def SFRoutput_plot(self, stream_name, nper=0, to_file=False):
    """TODO: Docs"""

    raise NotImplementedError("Tried to call unfinished method!")
    sfrout = SfrFile(os.path.join(self.data_folder, self.name + ".sfr.out"))
    sfr_df = sfrout.get_dataframe()
    sfr_df_time = sfr_df[sfr_df['time'] == nper]
    sfr_info = pd.DataFrame.from_records(self.model_data.boundaries.bc[stream_name]['bc_array'][0])
    sfr_info.loc[:, 'segment'] = sfr_info['iseg']
    sfr_df_time = pd.merge(sfr_df_time, sfr_info)
    sfr_df_time
# End SFRoutput_plot()


def plot_ts_obs_vs_sim(self, obs_name, obs_grp_name, obs_type):
    '''
    Plot a time series of the
    '''
    if obs_type == 'head':
        head_obs_ts = self.model_data.observations[obs_grp_name]['time_series']
        head_obs_ts[head_obs_ts['name'] == obs_name].plot(x='datetime', y='value')


def compareAllObs_metrics(self, to_file=False, head_name='head'):

    headobj = self.importHeads()
    times = headobj.get_times()

    scatterx = []
    scattery = []
    obs_sim_zone_all = []

    # The definition of obs_sim_zone looks like:
    # self.obs_sim_zone += [[obs, sim, zone, x, y]]

    for i in range(self.model_data.model_time.t['steps']):
        head = headobj.get_data(totim=times[i])
        self.CompareObserved(head_name, head, nper=i)
        obs_sim_zone_all += self.obs_sim_zone

    scatterx = np.array([h[0] for h in obs_sim_zone_all])
    scattery = np.array([h[1] for h in obs_sim_zone_all])

    ME = metric_me(scattery, scatterx)

    PBIAS = metric_pbias(scattery, scatterx)

    RMSE = metric_rmse(scattery, scatterx)

    SRMS = metric_srms(scattery, scatterx)

    if to_file:
        with open(os.path.join(self.data_folder, 'Head_Obs_Model_Measures.txt'), 'w') as f:
            f.write('ME PBIAS RMSE\n')
            f.write('{} {}% {} {}%'.format(ME, PBIAS, RMSE, SRMS))

    return ME, PBIAS, RMSE, SRMS
# End compareAllObs_metrics()


def compareAllObs(self, head_name):
    """TODO: Docs"""

    headobj = self.importHeads()
    times = headobj.get_times()

    obs_sim_zone_all = []

    # Could do it like this but need to doublecheck intent
    # for i, t in enumerate(times):
    #     head = headobj.get_data(totime=t)
    #
    #     # new method sets and returns `obs_sim_zone`
    #     obs_sim_zone_all += self.compare_observed('head', head, nper=i)
    #     # self.compare_observed('head', head, nper=i)
    #     # obs_sim_zone_all += self.obs_sim_zone
    # # End for

    # The definition of obs_sim_zone looks like:
    # self.obs_sim_zone += [[obs, sim, zone, x, y]]
    for i in range(self.model_data.model_time.t['steps']):
        head = headobj.get_data(totim=times[i])
        self.compare_observed(head_name, head, nper=i)
        obs_sim_zone_all += self.obs_sim_zone
    # End for

    scatterx = np.array([h[0] for h in obs_sim_zone_all])
    scattery = np.array([h[1] for h in obs_sim_zone_all])

    zoned_residuals = {}
    for i in range(1, 8):
        zoned_residuals[i] = [loc[0] - loc[1] for loc in obs_sim_zone_all if loc[2] == float(i)]

    residuals = [loc[0] - loc[1] for loc in obs_sim_zone_all]

    # First step is to set up the plot
    width = 20
    height = 5
    multiplier = 1.
    fig = plt.figure(figsize=(width * multiplier, height * multiplier))

    ax = fig.add_subplot(1, 3, 1)  # , aspect='equal')
    ax.set_title('Residuals')
    comp_zone_plots = {}
    colours = ['r', 'orangered', 'y', 'green', 'teal', 'blue', 'fuchsia']
    labels = ('qa', 'utb', 'utqa', 'utam', 'utaf', 'lta', 'bse')

    ax.hist(residuals, bins=20, alpha=0.5, color='black', histtype='step', label='all')
    for i in range(1, 8):
        comp_zone_plots[i] = ax.hist(zoned_residuals[i], bins=20, alpha=0.5,
                                     color=colours[i - 1], histtype='step', label=labels[i - 1])
    # End for

    plt.legend(loc='upper left', ncol=4, fontsize=11)
    ylim = ax.get_ylim()
    ax.set_ylim(ylim[0], ylim[1] + 0.25 * (ylim[1] - ylim[0]))

    ax = fig.add_subplot(1, 3, 2)
    ax.set_title('Sim vs Obs (%d points)' % (len(scatterx)))

    comp_zone_plots = {}
    for i in range(1, 8):
        scatterx2 = [loc[0] for loc in obs_sim_zone_all if loc[2] == float(i)]
        scattery2 = [loc[1] for loc in obs_sim_zone_all if loc[2] == float(i)]
        comp_zone_plots[i] = ax.scatter(scatterx2, scattery2, edgecolors=colours[
                                        i - 1], facecolors='none', alpha=0.5)

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

    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()

    ax.text(xmin + 0.45 * (xmax - xmin), ymin + 0.4 * (ymax - ymin),
            'Model Efficiency = %4.2f' % (metric_me(scattery, scatterx)))

    ax.text(xmin + 0.45 * (xmax - xmin), ymin + 0.3 * (ymax - ymin),
            'PBIAS = %4.2f%%' % (metric_pbias(scattery, scatterx)))

    ax.text(xmin + 0.45 * (xmax - xmin), ymin + 0.2 * (ymax - ymin),
            'RMSE = %4.2f' % (metric_rmse(scattery, scatterx)))

    #ax.plot(ax.get_ylim(), ax.get_ylim())

    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    new = (min(xlim[0], ylim[0]), max(xlim[1], ylim[1]))
#    new_upper = (min(xlim[0], ylim[0]) + unc, max(xlim[1], ylim[1]) + unc)
#    new_lower = (min(xlim[0], ylim[0]) - unc, max(xlim[1], ylim[1]) - unc)
    ax.plot(new, new, color='grey')
#    ax.plot(new_upper, new_upper, color='grey')
#    ax.plot(new_lower, new_lower, color='grey')
#    ax.fill_between(new, new_lower, new_upper, color='grey', alpha=0.3)
    ax.set_xlim(new)
    ax.set_ylim(new)

    ax = fig.add_subplot(1, 3, 3)
    ax.set_title('Residuals in space')

    modelmap = flopy.plot.ModelMap(model=self.mf)
    modelmap.plot_ibound()
    #modelmap.plot_bc('SFR', alpha=0.5, color='grey')
    x = np.array([h[3] for h in obs_sim_zone_all])
    y = np.array([h[4] for h in obs_sim_zone_all])
    zone = np.array([h[2] for h in obs_sim_zone_all])
    residuals = np.array([h[0] - h[1] for h in obs_sim_zone_all])

    colors_ = list(six.iteritems(colors.cnames))

    # Add the single letter colors.
    for name, rgb in six.iteritems(colors.ColorConverter.colors):
        colors_.append((name, colors.rgb2hex(rgb)))

    hex_ = [color[1] for color in colors_]
    nams = [color[0] for color in colors_]

    # Get the rgb equivalent.
    rgb_all = [colors.hex2color(color) for color in hex_]
    rgb_ref = []
    for col in colours:
        for index, nam in enumerate(nams):
            if col == nam:
                rgb_ref += [rgb_all[index]]

    rgba_colors = np.zeros((len(x), 4))

    # for red the first column needs to be one
    for i in range(1, 8):
        rgba_colors[:, 0][zone == i] = rgb_ref[i - 1][0]
        rgba_colors[:, 1][zone == i] = rgb_ref[i - 1][1]
        rgba_colors[:, 2][zone == i] = rgb_ref[i - 1][2]
    # End for

    # the fourth column needs to be your alphas
    rgba_colors[:, 3] = residuals / np.max(residuals)

    plt.scatter(x, y, c=residuals, alpha=0.5, edgecolors='none')
    start, end = ax.get_xlim()
    start = start // 1000 * 1000 + 1000
    end = end // 1000 * 1000 - 1000
    ax.xaxis.set_ticks(np.arange(start, end, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    self.plotRiverFromRiverSegData(ax, alpha=0.4)

    plt.colorbar()
    plt.show()
# End compareAllObs()


def compareAllObs2(self):

    # Set model output arrays to None to initialise
    head = None
    sfr_df = None
    stream_options = ['stage', 'depth', 'discharge']
    # Write observation to file
    for obs_set in self.model_data.observations.obs_group.keys():
        if self.model_data.observations.obs_group[obs_set]['real'] == False:
            continue
        # end if
        print("Processing {}".format(obs_set))

        obs_sim_zone_all = []

        obs_type = self.model_data.observations.obs_group[obs_set]['obs_type']
        # Import the required model outputs for processing
        if obs_type == 'head':
            # Check if model outputs have already been imported and if not import
            if not head:
                headobj = self.importHeads()
                head = headobj.get_alldata()
        elif obs_type in stream_options:
            try:
                sfr_df = self.sfr_df
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
            #sfr_chainage = self.model_data.observations.obs_group[obs_set]['locations']['Distance_Eppalock']
            for zone in obs_df['zone'].unique():
                if len(obs_df['zone'].unique()) == 1:
                    zone_txt = obs_set
                else:
                    zone_txt = obs_set + zone
                # End if
                obs_df_zone = obs_df[obs_df['zone'] == zone]
                for observation in obs_df_zone.index:
                    interval = int(obs_df_zone['interval'].loc[observation])
                    name = obs_df_zone['name'].loc[observation]
                    obs = obs_df_zone['value'].loc[observation]
                    #time = obs_df_zone['datetime'].loc[observation]
                    seg = sfr_location.loc[name]
                    #chainage = sfr_chainage.loc[name]
                    sfr = sfr_df
                    col_of_interest = obs_type
                    if obs_type == 'discharge':
                        col_of_interest = 'Qout'
                    sim_obs = sfr[(sfr['segment'] == seg) &
                                  (sfr['time'] == interval)][col_of_interest].tolist()[0]

                    obs_sim_zone_all += [[obs, sim_obs, seg]]  # , chainage, time]]

                # End for
            # End for
        # End if

        if obs_type == 'head':
            for zone in obs_df['zone'].unique():
                if len(obs_df['zone'].unique()) == 1:
                    zone_txt = obs_set #'head'
                else:
                    zone_txt = obs_set + zone
                # End if
                obs_df_zone = obs_df[obs_df['zone'] == zone]
                for observation in obs_df_zone.index:
                    interval = int(obs_df_zone['interval'].loc[observation])
                    name = obs_df_zone['name'].loc[observation]
                    obs = obs_df_zone['value'].loc[observation]
                    (x_cell, y_cell) = self.model_data.mesh2centroid2Dindex[
                        (sim_map_dict[name][1], sim_map_dict[name][2])]
                    (lay, row, col) = [sim_map_dict[name][0],
                                       sim_map_dict[name][1], sim_map_dict[name][2]]

                    sim_heads = [head[interval][lay][row][col]]

                    sim_head = np.mean(sim_heads)  # ???
                    obs_sim_zone_all += [[obs, sim_head, zone]]

                # End for
            # End for
        # End if

        plot_obs_vs_sim(obs_set, obs_sim_zone_all, unc=2)


# End compareAllObs()

def plot_bc(self, modelmap, name):  # , **kwargs):
    '''
    Plot boundary condition in active cells of passed modelmap object.
    Includes rough test for overcoming non-existence of bc
    '''
    try:
        modelmap.plot_bc(name)  # , **kwargs)
    except:
        print("Could not find '{}' package in model".format(name))


def plot_bcs(self, modelmap, bcs):  # , **kwargs):
    for bc in bcs:
        self.plot_bc(modelmap, bc)  # , **kwargs)


def plot_bc_by_layer(self, bc_name, plot_rows=2):
    
    num_plots = self.nlay
    plot_cols = num_plots / plot_rows + 1
    
    fig = plt.figure()
    for plot in range(num_plots):
        fig.add_subplot(plot_rows, plot_cols, plot + 1, aspect='equal')
        modelmap = flopy.plot.ModelMap(model=self.mf, layer=plot)
        modelmap.plot_ibound()
        self.plot_bc(modelmap, bc_name)
        
        
def viewHeadsByZone(self, nper='all', head_name='head'):
    """
    :param nper: Default value = 'all')
    """

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
        self.compare_observed(head_name, head_orig, nper=nper)
    # End if

    # First step is to set up the plot
    width = 20
    height = 10
    multiplier = 1.0
    fig = plt.figure(figsize=(width * multiplier, height * multiplier))

    vmin = 0
    vmax = 200

    ax = fig.add_subplot(2, 4, 1, aspect='equal')

    ax.set_title('ibound and bc')
    # Next we create an instance of the ModelMap class
    modelmap = flopy.plot.ModelMap(model=self.mf)
    modelmap.plot_ibound()

    bcs = ['RIV', 'SFR', 'GHB']
    self.plot_bcs(modelmap, bcs)

    ax.axes.xaxis.set_ticklabels([])

    ax = fig.add_subplot(2, 4, 2, aspect='equal')
    ax.set_title('Coonambidgal')
    modelmap = flopy.plot.ModelMap(model=self.mf)
    max_head = np.amax(head)
    min_head = np.amin(head)

    modelmap_plot_values = {
        'masked_values': [-999.98999023, max_head, min_head, np.nan],
        'alpha': 0.5,
        'vmin': vmin,
        'vmax': vmax
    }

    array = modelmap.plot_array(head[0], **modelmap_plot_values)

    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])
    cbar_ax2 = fig.add_axes([0.43, 0.525, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax2)

    ax = fig.add_subplot(2, 4, 3, aspect='equal')
    ax.set_title('Shepparton')
    modelmap = flopy.plot.ModelMap(model=self.mf)
    array = modelmap.plot_array(head[2], **modelmap_plot_values)
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])
    cbar_ax1 = fig.add_axes([0.67, 0.525, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax1)

    ax = fig.add_subplot(2, 4, 5, aspect='equal')
    ax.set_title('Calivil')
    modelmap = flopy.plot.ModelMap(model=self.mf)
    array = modelmap.plot_array(head[4], **modelmap_plot_values)
    start, end = ax.get_xlim()
    start = start // 1000 * 1000 + 1000
    end = end // 1000 * 1000 - 1000
    ax.xaxis.set_ticks(np.arange(start, end, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    cbar_ax3 = fig.add_axes([0.19, 0.055, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax3)

    ax = fig.add_subplot(2, 4, 6, aspect='equal')
    ax.set_title('Renmark')
    modelmap = flopy.plot.ModelMap(model=self.mf)
    array = modelmap.plot_array(head[5], **modelmap_plot_values)
    ax.yaxis.set_ticklabels([])
    ax.xaxis.set_ticks(np.arange(start, end, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    cbar_ax4 = fig.add_axes([0.43, 0.055, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax4)

    ax = fig.add_subplot(2, 4, 7, aspect='equal')
    ax.set_title('Basement')
    modelmap = flopy.plot.ModelMap(model=self.mf)
    array = modelmap.plot_array(head[6], **modelmap_plot_values)
    ax.yaxis.set_ticklabels([])
    ax.xaxis.set_ticks(np.arange(start, end, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    cbar_ax5 = fig.add_axes([0.67, 0.055, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax5)

    ax.plot(ax.get_ylim(), ax.get_ylim())
    fig.subplots_adjust(left=0.01, right=0.95, bottom=0.05, top=0.95, wspace=0.1, hspace=0.12)
    plt.show()

# End viewHeadsByZone()


def viewHeadsByZone2(self, iter_num, nper='all', head_name='head'):
    """
    :param iter_num: param nper:  (Default value = 'all')

    :param nper:  (Default value = 'all')
    """

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
    # End if

    if nper == 'all':
        scatterx = []
        scattery = []
        obs_sim_zone_all = []
        for i in range(len(times)):
            self.compare_observed(head_name, head_orig, nper=i)
            scatterx += [h[0] for h in self.obs_sim_zone]
            scattery += [h[1] for h in self.obs_sim_zone]
            obs_sim_zone_all += self.obs_sim_zone
        self.obs_sim_zone = obs_sim_zone_all
    else:
        self.compare_observed(head_name, head_orig, nper=nper)
        scatterx = [h[0] for h in self.obs_sim_zone]
        scattery = [h[1] for h in self.obs_sim_zone]
    # End if

    # First step is to set up the plot
    width = 20
    height = 10
    multiplier = 1.0
    fig = plt.figure(figsize=(width * multiplier, height * multiplier))

    vmin = 0
    vmax = 200

    ax = fig.add_subplot(2, 4, 1, aspect='equal')

    ax.set_title('ibound and bc')
    # Next we create an instance of the ModelMap class
    modelmap = flopy.plot.ModelMap(model=self.mf)
    modelmap.plot_ibound()

    bcs = ['RIV', 'WEL', 'SFR', 'GHB']
    self.plot_bcs(modelmap, bcs)

    start, end = ax.get_xlim()
    start = start // 1000 * 1000 + 1000
    end = end // 1000 * 1000 - 1000
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))

    ax = fig.add_subplot(2, 4, 2, aspect='equal')
    ax.set_title('Coonambidgal')
    modelmap = flopy.plot.ModelMap(model=self.mf)  # , sr=self.mf.dis.sr, dis=self.mf.dis)
    max_head = np.amax(head)
    min_head = np.amin(head)

    array = modelmap.plot_array(
        head[0], masked_values=[-999.98999023, max_head, min_head], alpha=0.5, vmin=vmin, vmax=vmax)

    ax.xaxis.set_ticks(np.arange(start, end, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    cbar_ax2 = fig.add_axes([0.43, 0.525, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax2)

    scatterx2 = [loc[3] for loc in self.obs_sim_zone if loc[2] == 1.0]
    scattery2 = [loc[4] for loc in self.obs_sim_zone if loc[2] == 1.0]
    ax.scatter(scatterx2, scattery2, c=[loc[0] for loc in self.obs_sim_zone if loc[
               2] == 1.0], alpha=0.8, vmin=vmin, vmax=vmax)

    ax.text(2.7e5, 6030000, 'Observations: %d' % (len(scatterx2)))

    ax = fig.add_subplot(2, 4, 3, aspect='equal')
    ax.set_title('Shepparton')
    modelmap = flopy.plot.ModelMap(model=self.mf)  # , sr=self.mf.dis.sr, dis=self.mf.dis)
    array = modelmap.plot_array(
        head[2], masked_values=[-999.98999023, max_head, min_head, np.nan], alpha=0.5, vmin=vmin, vmax=vmax)

    ax.xaxis.set_ticks(np.arange(start, end, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    cbar_ax1 = fig.add_axes([0.67, 0.525, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax1)

    scatterx2 = [loc[3] for loc in self.obs_sim_zone if loc[2] == 3.0]
    scattery2 = [loc[4] for loc in self.obs_sim_zone if loc[2] == 3.0]
    ax.scatter(scatterx2, scattery2, c=[loc[0] for loc in self.obs_sim_zone if loc[
               2] == 3.0], alpha=0.8, vmin=vmin, vmax=vmax)
    ax.text(2.7e5, 6030000, 'Observations: %d' % (len(scatterx2)))

    ax = fig.add_subplot(2, 4, 4)
    ax.set_title('Residuals')
    ax.hist([loc[0] - loc[1] for loc in self.obs_sim_zone], bins=20, alpha=0.5)

    ax = fig.add_subplot(2, 4, 5, aspect='equal')
    ax.set_title('Calivil')
    modelmap = flopy.plot.ModelMap(model=self.mf)
    array = modelmap.plot_array(
        head[4], masked_values=[-999.98999023, max_head, min_head, np.nan], alpha=0.5, vmin=vmin, vmax=vmax)
    ax.xaxis.set_ticks(np.arange(start, end, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    cbar_ax3 = fig.add_axes([0.19, 0.055, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax3)

    scatterx2 = [loc[3] for loc in self.obs_sim_zone if loc[2] == 5.0]
    scattery2 = [loc[4] for loc in self.obs_sim_zone if loc[2] == 5.0]
    ax.scatter(scatterx2, scattery2, c=[loc[0] for loc in self.obs_sim_zone if loc[
               2] == 5.0], alpha=0.8, vmin=vmin, vmax=vmax)
    ax.text(2.7e5, 6030000, 'Observations: %d' % (len(scatterx2)))

    ax = fig.add_subplot(2, 4, 6, aspect='equal')
    ax.set_title('Renmark')
    modelmap = flopy.plot.ModelMap(model=self.mf)
    array = modelmap.plot_array(
        head[5], masked_values=[-999.98999023, max_head, min_head, np.nan], alpha=0.5, vmin=vmin, vmax=vmax)

    ax.xaxis.set_ticks(np.arange(start, end, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    cbar_ax4 = fig.add_axes([0.43, 0.055, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax4)

    scatterx2 = [loc[3] for loc in self.obs_sim_zone if loc[2] == 6.0]
    scattery2 = [loc[4] for loc in self.obs_sim_zone if loc[2] == 6.0]
    ax.scatter(scatterx2, scattery2, c=[loc[0] for loc in self.obs_sim_zone if loc[
               2] == 6.0], alpha=0.8, vmin=vmin, vmax=vmax)
    ax.text(2.7e5, 6030000, 'Observations: %d' % (len(scatterx2)))

    ax = fig.add_subplot(2, 4, 7, aspect='equal')
    ax.set_title('Basement')
    modelmap = flopy.plot.ModelMap(model=self.mf)
    array = modelmap.plot_array(
        head[6], masked_values=[-999.98999023, max_head, min_head, np.nan], alpha=0.5, vmin=vmin, vmax=vmax)

    ax.xaxis.set_ticks(np.arange(start, end, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    cbar_ax5 = fig.add_axes([0.67, 0.055, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax5)

    scatterx2 = [loc[3] for loc in self.obs_sim_zone if loc[2] == 7.0]
    scattery2 = [loc[4] for loc in self.obs_sim_zone if loc[2] == 7.0]
    ax.scatter(scatterx2, scattery2, c=[loc[0] for loc in self.obs_sim_zone if loc[
               2] == 7.0], alpha=0.8, vmin=vmin, vmax=vmax)
    ax.text(2.7e5, 6030000, 'Observations: %d' % (len(scatterx2)))

    ax = fig.add_subplot(2, 4, 8, aspect=0.9)
    ax.set_title('Sim vs Obs (%d points)' % (len(scatterx)))

    comp_zone_plots = {}
    colors = ['b', 'c', 'y', 'm', 'r', 'green', 'orange']
    for i in range(1, 8):
        scatterx2 = [loc[0] for loc in self.obs_sim_zone if loc[2] == float(i)]
        scattery2 = [loc[1] for loc in self.obs_sim_zone if loc[2] == float(i)]
        comp_zone_plots[i] = ax.scatter(scatterx2, scattery2, edgecolors=colors[
                                        i - 1], facecolors='none', alpha=0.5)

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

    ax.text(150, 75, 'Model Efficiency = %4.2f' % (metric_me(scattery, scatterx)))

    ax.text(150, 40, 'PBIAS = %4.2f%%' % (metric_pbias(scattery, scatterx)))

    ax.text(150, 20, 'RMSE = %4.2f' % (metric_rmse(scattery, scatterx)))

    ax.plot(ax.get_ylim(), ax.get_ylim())

    fig.subplots_adjust(left=0.01, right=0.95, bottom=0.05, top=0.95, wspace=0.1, hspace=0.12)
    plt.savefig('run_viewheads_{}.png'.format(iter_num), bbox_inches='tight')

# End viewHeadsByZone2()


def viewHeadsByZone3(self, iter_num, nper='all', head_name='head'):
    """
    :param iter_num: param nper:  (Default value = 'all')

    :param nper:  (Default value = 'all')
    """

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
    # End if

    if nper == 'all':
        scatterx = []
        scattery = []
        obs_sim_zone_all = []
        for i in range(len(times)):
            self.compare_observed(head_name, head_orig, nper=i)
            scatterx += [h[0] for h in self.obs_sim_zone]
            scattery += [h[1] for h in self.obs_sim_zone]
            obs_sim_zone_all += self.obs_sim_zone
        self.obs_sim_zone = obs_sim_zone_all
    else:
        self.compare_observed(head_name, head_orig, nper=nper)
        scatterx = [h[0] for h in self.obs_sim_zone]
        scattery = [h[1] for h in self.obs_sim_zone]
    # End if

    # First step is to set up the plot
    width = 20
    height = 10
    multiplier = 1.0
    fig = plt.figure(figsize=(width * multiplier, height * multiplier))

    vmin = 0
    vmax = 200

    ax = fig.add_subplot(2, 4, 1, aspect='equal')

    ax.set_title('ibound and bc')
    # Next we create an instance of the ModelMap class
    modelmap = flopy.plot.ModelMap(model=self.mf)
    modelmap.plot_ibound()

    bcs = ['RIV', 'WEL', 'SFR', 'GHB']
    self.plot_bcs(modelmap, bcs)

    start, end = ax.get_xlim()
    start = start // 1000 * 1000 + 1000
    end = end // 1000 * 1000 - 1000
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))

    ax = fig.add_subplot(2, 4, 2, aspect='equal')
    ax.set_title('Coonambidgal')
    modelmap = flopy.plot.ModelMap(model=self.mf)  # , sr=self.mf.dis.sr, dis=self.mf.dis)
    max_head = np.amax(head)
    min_head = np.amin(head)

    array = modelmap.plot_array(
        head[0], masked_values=[-999.98999023, max_head, min_head], alpha=0.5, vmin=vmin, vmax=vmax)

    ax.xaxis.set_ticks(np.arange(start, end, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    cbar_ax2 = fig.add_axes([0.43, 0.525, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax2)

    scatterx2 = [loc[3] for loc in self.obs_sim_zone if loc[2] == 1.0]
    scattery2 = [loc[4] for loc in self.obs_sim_zone if loc[2] == 1.0]
    ax.scatter(scatterx2, scattery2, c=[loc[0] - loc[1] for loc in self.obs_sim_zone if loc[
               2] == 1.0], alpha=0.8, vmin=vmin, vmax=vmax)

    ax.text(2.7e5, 6030000, 'Observations: %d' % (len(scatterx2)))

    ax = fig.add_subplot(2, 4, 3, aspect='equal')
    ax.set_title('Shepparton')
    modelmap = flopy.plot.ModelMap(model=self.mf)  # , sr=self.mf.dis.sr, dis=self.mf.dis)
    array = modelmap.plot_array(
        head[2], masked_values=[-999.98999023, max_head, min_head, np.nan], alpha=0.5, vmin=vmin, vmax=vmax)

    ax.xaxis.set_ticks(np.arange(start, end, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    cbar_ax1 = fig.add_axes([0.67, 0.525, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax1)

    scatterx2 = [loc[3] for loc in self.obs_sim_zone if loc[2] == 3.0]
    scattery2 = [loc[4] for loc in self.obs_sim_zone if loc[2] == 3.0]
    ax.scatter(scatterx2, scattery2, c=[loc[0] - loc[1] for loc in self.obs_sim_zone if loc[
               2] == 3.0], alpha=0.8, vmin=vmin, vmax=vmax)
    ax.text(2.7e5, 6030000, 'Observations: %d' % (len(scatterx2)))

    ax = fig.add_subplot(2, 4, 4)
    ax.set_title('Residuals')
    ax.hist([loc[0] - loc[1] for loc in self.obs_sim_zone], bins=20, alpha=0.5)

    ax = fig.add_subplot(2, 4, 5, aspect='equal')
    ax.set_title('Calivil')
    modelmap = flopy.plot.ModelMap(model=self.mf)
    array = modelmap.plot_array(
        head[4], masked_values=[-999.98999023, max_head, min_head, np.nan], alpha=0.5, vmin=vmin, vmax=vmax)
    ax.xaxis.set_ticks(np.arange(start, end, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    cbar_ax3 = fig.add_axes([0.19, 0.055, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax3)

    scatterx2 = [loc[3] for loc in self.obs_sim_zone if loc[2] == 5.0]
    scattery2 = [loc[4] for loc in self.obs_sim_zone if loc[2] == 5.0]
    ax.scatter(scatterx2, scattery2, c=[loc[0] - loc[1] for loc in self.obs_sim_zone if loc[
               2] == 5.0], alpha=0.8, vmin=vmin, vmax=vmax, cmap='viridis')
    ax.text(2.7e5, 6030000, 'Observations: %d' % (len(scatterx2)))

    ax = fig.add_subplot(2, 4, 6, aspect='equal')
    ax.set_title('Renmark')
    modelmap = flopy.plot.ModelMap(model=self.mf)
    array = modelmap.plot_array(
        head[5], masked_values=[-999.98999023, max_head, min_head, np.nan], alpha=0.5, vmin=vmin, vmax=vmax)

    ax.xaxis.set_ticks(np.arange(start, end, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    cbar_ax4 = fig.add_axes([0.43, 0.055, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax4)

    scatterx2 = [loc[3] for loc in self.obs_sim_zone if loc[2] == 6.0]
    scattery2 = [loc[4] for loc in self.obs_sim_zone if loc[2] == 6.0]
    ax.scatter(scatterx2, scattery2, c=[loc[0] - loc[1] for loc in self.obs_sim_zone if loc[
               2] == 6.0], alpha=0.8, vmin=vmin, vmax=vmax)
    ax.text(2.7e5, 6030000, 'Observations: %d' % (len(scatterx2)))

    ax = fig.add_subplot(2, 4, 7, aspect='equal')
    ax.set_title('Basement')
    modelmap = flopy.plot.ModelMap(model=self.mf)
    array = modelmap.plot_array(
        head[6], masked_values=[-999.98999023, max_head, min_head, np.nan], alpha=0.5, vmin=vmin, vmax=vmax)

    ax.xaxis.set_ticks(np.arange(start, end, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    cbar_ax5 = fig.add_axes([0.67, 0.055, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax5)

    scatterx2 = [loc[3] for loc in self.obs_sim_zone if loc[2] == 7.0]
    scattery2 = [loc[4] for loc in self.obs_sim_zone if loc[2] == 7.0]
    ax.scatter(scatterx2, scattery2, c=[loc[0] - loc[1] for loc in self.obs_sim_zone if loc[
               2] == 7.0], alpha=0.8, vmin=vmin, vmax=vmax)
    ax.text(2.7e5, 6030000, 'Observations: %d' % (len(scatterx2)))

    ax = fig.add_subplot(2, 4, 8, aspect=0.9)
    ax.set_title('Sim vs Obs (%d points)' % (len(scatterx)))

    comp_zone_plots = {}
    colors = ['b', 'c', 'y', 'm', 'r', 'green', 'orange']
    for i in range(1, 8):
        scatterx2 = [loc[0] for loc in self.obs_sim_zone if loc[2] == float(i)]
        scattery2 = [loc[1] for loc in self.obs_sim_zone if loc[2] == float(i)]
        comp_zone_plots[i] = ax.scatter(scatterx2, scattery2, edgecolors=colors[
                                        i - 1], facecolors='none', alpha=0.5)

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

    ax.text(150, 75, 'Model Efficiency = %4.2f' % (metric_me(scattery, scatterx)))

    ax.text(150, 60, 'PBIAS = %4.2f%%' % (metric_pbias(scattery, scatterx)))

    ax.text(150, 45, 'RMSE = %4.2f' % (metric_rmse(scattery, scatterx)))

    ax.plot(ax.get_ylim(), ax.get_ylim())
    fig.subplots_adjust(left=0.01, right=0.95, bottom=0.05, top=0.95, wspace=0.1, hspace=0.12)
    plt.savefig('run_viewheads_{}.png'.format(iter_num), bbox_inches='tight')

# End viewHeadsByZone3()


def viewHeads(self):
    """TODO: Docs"""

    # Create the headfile object
    headobj = self.importHeads()
    # cbbobj = self.import_cbb()
    times = headobj.get_times()
    head = headobj.get_data(totim=times[-1])

    # First step is to set up the plot
    fig = plt.figure(figsize=(20, 10))
    ax = fig.add_subplot(2, 4, 1, aspect='equal')
    ax.set_title('ibound and bc')
    # Next we create an instance of the ModelMap class
    modelmap = flopy.plot.ModelMap(model=self.mf)
    modelmap.plot_ibound()
    self.plotRiverFromRiverSegData(ax)
    bcs = ['RIV', 'SFR']
    self.plot_bcs(modelmap, bcs)
    ax.axes.xaxis.set_ticklabels([])
    vmin = np.round(np.amin(head[head != -999.99]))
    vmax = np.amax(head)
    cmap = 'jet'

    ax = fig.add_subplot(2, 4, 2, aspect='equal')
    ax.set_title('Heads layer 1')
    modelmap = flopy.plot.ModelMap(model=self.mf)
    modelmap.plot_ibound()
    max_head = np.amax(head)
    min_head = np.amin(head)

    levels = np.arange(vmin, vmax, 10)
    array = modelmap.contour_array(head[0], masked_values=[-999.98999023, max_head,
                                                           min_head], alpha=0.9, vmin=vmin, vmax=vmax, cmap=cmap, levels=levels)
    plt.clabel(array, inline=1, fontsize=10)
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])
    self.plotRiverFromRiverSegData(ax)
    cbar_ax2 = fig.add_axes([0.43, 0.525, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax2)

    ax = fig.add_subplot(2, 4, 3, aspect='equal')
    ax.set_title('Heads layer 2')
    modelmap = flopy.plot.ModelMap(model=self.mf)
    modelmap.plot_ibound()
    array = modelmap.contour_array(head[1], masked_values=[-999.98999023, max_head,
                                                           min_head], alpha=0.9, vmin=vmin, vmax=vmax, cmap=cmap, levels=levels)
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])
    self.plotRiverFromRiverSegData(ax)
    cbar_ax1 = fig.add_axes([0.67, 0.525, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax1)

    ax = fig.add_subplot(2, 4, 4, aspect='equal')
    ax.set_title('Heads layer 3')
    modelmap = flopy.plot.ModelMap(model=self.mf)
    modelmap.plot_ibound()
    array = modelmap.contour_array(head[2], masked_values=[-999.98999023, max_head,
                                                           min_head], alpha=0.9, vmin=vmin, vmax=vmax, cmap=cmap, levels=levels)
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])
    self.plotRiverFromRiverSegData(ax)
    cbar_ax5 = fig.add_axes([0.91, 0.525, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax5)

    ax = fig.add_subplot(2, 4, 5, aspect='equal')
    ax.set_title('Heads layer 4')
    modelmap = flopy.plot.ModelMap(model=self.mf)
    modelmap.plot_ibound()
    array = modelmap.contour_array(head[3], masked_values=[-999.98999023, max_head,
                                                           min_head], alpha=0.9, vmin=vmin, vmax=vmax, cmap=cmap, levels=levels)
    start, end = ax.get_xlim()
    start = start // 1000 * 1000 + 1000
    end = end // 1000 * 1000 - 1000
    ax.xaxis.set_ticks(np.arange(start, end, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    self.plotRiverFromRiverSegData(ax)

    cbar_ax3 = fig.add_axes([0.19, 0.055, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax3)

    ax = fig.add_subplot(2, 4, 6, aspect='equal')
    ax.set_title('Heads layer 5')
    modelmap = flopy.plot.ModelMap(model=self.mf)
    modelmap.plot_ibound()
    array = modelmap.contour_array(head[4], masked_values=[-999.98999023, max_head,
                                                           min_head], alpha=0.9, vmin=vmin, vmax=vmax, cmap=cmap, levels=levels)
    ax.yaxis.set_ticklabels([])
    ax.xaxis.set_ticks(np.arange(start, end, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    self.plotRiverFromRiverSegData(ax)

    cbar_ax4 = fig.add_axes([0.43, 0.055, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax4)

    ax = fig.add_subplot(2, 4, 7, aspect='equal')
    ax.set_title('Heads layer 6')
    modelmap = flopy.plot.ModelMap(model=self.mf)
    modelmap.plot_ibound()
    array = modelmap.contour_array(head[5], masked_values=[-999.98999023, max_head,
                                                           min_head], alpha=0.9, vmin=vmin, vmax=vmax, cmap=cmap, levels=levels)
    self.plotRiverFromRiverSegData(ax)
    ax.yaxis.set_ticklabels([])
    ax.xaxis.set_ticks(np.arange(start, end, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))

    cbar_ax5 = fig.add_axes([0.67, 0.055, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax5)
    fig.subplots_adjust(left=0.01, right=0.95, bottom=0.05, top=0.95, wspace=0.1, hspace=0.12)
    plt.show()

# End viewHeads()


def viewHeadLayer(self, layer=0, figsize=(20, 10)):
    """
    :param layer: (Default value = 0)

    :param figsize: (Default value = (20, 10))
    """
    # Create the headfile object
    headobj = self.importHeads()
    #cbbobj = self.importCbb()
    times = headobj.get_times()
    head = headobj.get_data(totim=times[-1])

    #frf = cbbobj.get_data(text='FLOW RIGHT FACE')[0]
    #fff = cbbobj.get_data(text='FLOW FRONT FACE')[0]

    # First step is to set up the plot
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(1, 2, 1, aspect='equal')
    ax.set_title('ibound and bc')
    # Next we create an instance of the ModelMap class
    modelmap = flopy.plot.ModelMap(model=self.mf)  # , sr=self.mf.dis.sr, dis=self.mf.dis)
    modelmap.plot_ibound()
    bcs = ['RIV']
    self.plot_bcs(modelmap, bcs)
    ax.axes.xaxis.set_ticklabels([])
    vmin = np.round(np.amin(head[head != -999.99]))
    vmax = 150
    cmap = 'jet'

    ax = fig.add_subplot(1, 2, 2, aspect='equal')
    ax.set_title('Heads layer {}'.format(layer + 1))
    modelmap = flopy.plot.ModelMap(model=self.mf)
    modelmap.plot_ibound()
    max_head = np.amax(head)
    min_head = np.amin(head)

    levels = np.arange(vmin, vmax, 1)

    modelmap.plot_bc('RIV', alpha=0.3, plotAll=True)
    array = modelmap.contour_array(head[0], masked_values=[-999.98999023, max_head,
                                                           min_head], alpha=0.9, vmin=vmin, vmax=vmax, cmap=cmap, levels=levels)
    plt.clabel(array, inline=True, fontsize=10)
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])
    cbar_ax2 = fig.add_axes([0.43, 0.525, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax2)

    fig.subplots_adjust(left=0.01, right=0.95, bottom=0.05, top=0.95, wspace=0.1, hspace=0.12)

    plt.show()

# End viewHeadLayer()


def viewHeads2(self):
    """TODO: Docs"""

    # Create the headfile object
    headobj = bf.HeadFile(os.path.join(self.data_folder, self.name + '.hds'))
    cbbobj = bf.CellBudgetFile(os.path.join(self.data_folder, self.name + '.cbc'))

    times = headobj.get_times()
    head = headobj.get_data(totim=times[0])

    water_balance_components = cbbobj.textlist
    water_balance = {}
    water_balance_summary = {}
    water_bal_summed_titles = []
    water_bal_summed_values = []
    for component in water_balance_components:
        component_stripped = component.lstrip().rstrip()
        if 'FLOW' in component_stripped:
            continue
        water_balance[component_stripped] = cbbobj.get_data(text=component, full3D=True)[-1]
        if np.any(water_balance[component_stripped][0] > 0):
            water_balance_summary[component_stripped + '_pos'] = np.sum(
                water_balance[component_stripped][0][water_balance[component_stripped][0] > 0])
        else:
            water_balance_summary[component_stripped + '_pos'] = 0.0
        if np.any(water_balance[component_stripped][0] < 0):
            water_balance_summary[component_stripped + '_neg'] = np.sum(
                water_balance[component_stripped][0][water_balance[component_stripped][0] < 0])
        else:
            water_balance_summary[component_stripped + '_neg'] = 0.0
        # End if

        water_bal_summed_titles += component_stripped + '_pos', component_stripped + '_neg'
        water_bal_summed_values += water_balance_summary[
            component_stripped + '_pos'], water_balance_summary[component_stripped + '_neg']
    # End for

    water_bal_summed_titles += ['Error']
    Error = sum(water_bal_summed_values)
    water_bal_summed_values += [Error]
    wat_bal_df = pd.DataFrame(water_bal_summed_values, water_bal_summed_titles)
    wat_bal_df.columns = ['Flux m^3/d']
    wat_bal_df = wat_bal_df[wat_bal_df['Flux m^3/d'] != 0.]

    # First step is to set up the plot
    fig = plt.figure(figsize=(20, 10))
    ax = fig.add_subplot(2, 4, 1, aspect='equal')
    ax.set_title('ibound and bc')
    # Next we create an instance of the ModelMap class
    modelmap = flopy.plot.ModelMap(model=self.mf)  # , sr=self.mf.dis.sr, dis=self.mf.dis)
    modelmap.plot_ibound()
    bcs = ['RIV', 'SFR', 'GHB']
    self.plot_bcs(modelmap, bcs)
    self.plotRiverFromRiverSegData(ax)
    ax.axes.xaxis.set_ticklabels([])

    ax = fig.add_subplot(2, 4, 2, aspect='equal')
    ax.set_title('Heads')
    modelmap = flopy.plot.ModelMap(model=self.mf)  # , sr=self.mf.dis.sr, dis=self.mf.dis)
    max_head = np.amax(head)
    min_head = np.amin(head)
    array = modelmap.plot_array(
        head, masked_values=[-999.98999023, max_head, min_head], alpha=0.5)
    self.plotRiverFromRiverSegData(ax)
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])
    cbar_ax2 = fig.add_axes([0.43, 0.525, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax2)

    ax = fig.add_subplot(2, 4, 3, aspect='equal')
    ax.set_title('Murray exchange')
    modelmap.plot_ibound()
    modelmap = flopy.plot.ModelMap(model=self.mf)  # , sr=self.mf.dis.sr, dis=self.mf.dis)
    river_flux = modelmap.plot_array(water_balance['RIVER LEAKAGE'].sum(axis=0), alpha=0.5)
    self.plotRiverFromRiverSegData(ax)
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])

    cbar_ax1 = fig.add_axes([0.67, 0.525, 0.01, 0.42])
    fig.colorbar(river_flux, cax=cbar_ax1)

    # Setup params to get water balance aspect ratio looking nice
    aspect = float(
        12.5715 / ((wat_bal_df.max()[0] - wat_bal_df.min()[0]) / float(wat_bal_df.shape[1])))

    ax = fig.add_subplot(2, 4, 4, aspect=aspect)
    modelmap.plot_ibound()
    ax.set_title('Campaspe exchange')
    modelmap = flopy.plot.ModelMap(model=self.mf)  # , sr=self.mf.dis.sr, dis=self.mf.dis)
    try:
        river_flux = modelmap.plot_array(water_balance['STREAM LEAKAGE'].sum(axis=0), alpha=0.5)
    except:
        river_flux = modelmap.plot_array(water_balance['RIVER LEAKAGE'].sum(axis=0), alpha=0.5)
    self.plotRiverFromRiverSegData(ax)
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])

    cbar_ax1 = fig.add_axes([0.92, 0.525, 0.01, 0.42])
    fig.colorbar(river_flux, cax=cbar_ax1)

    # Setup params to get water balance aspect ratio looking nice
    aspect = float(
        12.5715 / ((wat_bal_df.max()[0] - wat_bal_df.min()[0]) / float(wat_bal_df.shape[1])))

    ax = fig.add_subplot(2, 4, 8, aspect=aspect)
    modelmap.plot_ibound()
    ax.set_title('GW under Murray')
    modelmap = flopy.plot.ModelMap(model=self.mf)  # , sr=self.mf.dis.sr, dis=self.mf.dis)
    modelmap.plot_ibound()
    river_flux = modelmap.plot_array(water_balance['HEAD DEP BOUNDS'].sum(axis=0), alpha=0.5)
    self.plotRiverFromRiverSegData(ax)
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])

    cbar_ax1 = fig.add_axes([0.92, 0.055, 0.01, 0.42])
    fig.colorbar(river_flux, cax=cbar_ax1)

    # Setup params to get water balance aspect ratio looking nice
    aspect = float(
        12.5715 / ((wat_bal_df.max()[0] - wat_bal_df.min()[0]) / float(wat_bal_df.shape[1])))

    ax = fig.add_subplot(2, 4, 5, aspect='equal')
    ax.set_title('Rainfall recharge')
    modelmap = flopy.plot.ModelMap(model=self.mf)  # , sr=self.mf.dis.sr, dis=self.mf.dis)
    modelmap.plot_ibound()
    recharge = modelmap.plot_array(water_balance['RECHARGE'][0], masked_values=[0.], alpha=0.5)
    self.plotRiverFromRiverSegData(ax)

    start, end = ax.get_xlim()
    start = start // 1000 * 1000 + 1000
    end = end // 1000 * 1000 - 1000
    ax.xaxis.set_ticks(np.arange(start, end, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))

    cbar_ax3 = fig.add_axes([0.19, 0.055, 0.01, 0.42])
    fig.colorbar(recharge, cax=cbar_ax3)

    ax = fig.add_subplot(2, 4, 6, aspect='equal')
    ax.set_title('Surface elevation')
    modelmap = flopy.plot.ModelMap(model=self.mf)  # , sr=self.mf.dis.sr, dis=self.mf.dis)
    modelmap.plot_ibound()
    elev = self.mf.dis.top.array
    elev = np.ma.masked_where(self.model_data.model_mesh3D[1][0] == -1, elev)
    elevation = modelmap.plot_array(elev, alpha=0.5)
    self.plotRiverFromRiverSegData(ax)
    ax.yaxis.set_ticklabels([])
    ax.xaxis.set_ticks(np.arange(start, end, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))

    cbar_ax4 = fig.add_axes([0.43, 0.055, 0.01, 0.42])
    fig.colorbar(elevation, cax=cbar_ax4)
    ax = fig.add_subplot(2, 4, 7, aspect='equal')
    ax.set_title('Pressure Head')
    modelmap = flopy.plot.ModelMap(model=self.mf)  # , sr=self.mf.dis.sr, dis=self.mf.dis)
    modelmap.plot_ibound()
    pressure = head[0] - elev
    storage = modelmap.plot_array(
        np.ma.masked_where(self.model_data.model_mesh3D[1][0] == -1, pressure), masked_values=[x for x in np.unique(pressure) if x > 0.], alpha=0.5)  # , vmin=-50, vmax=50)
    self.plotRiverFromRiverSegData(ax)
    ax.yaxis.set_ticklabels([])
    ax.xaxis.set_ticks(np.arange(start, end, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))

    cbar_ax5 = fig.add_axes([0.67, 0.055, 0.01, 0.42])
    fig.colorbar(storage, cax=cbar_ax5)

    fig.subplots_adjust(left=0.01, right=0.95, bottom=0.05, top=0.95,
                        wspace=0.1, hspace=0.12)
    plt.show()

    cbbobj.close()
    headobj.close()

# End viewHeads2()


def viewGHB(self):
    """TODO: Docs"""

    # Create the headfile object
    headobj = bf.HeadFile(os.path.join(self.data_folder, self.name + '.hds'))
    cbbobj = bf.CellBudgetFile(os.path.join(self.data_folder, self.name + '.cbc'))

    water_balance_components = cbbobj.textlist
    water_balance = {}
    water_balance_summary = {}
    water_bal_summed_titles = []
    water_bal_summed_values = []

    for component in water_balance_components:
        component_stripped = component.lstrip().rstrip()
        if 'FLOW' in component_stripped:
            continue
        # End if
        water_balance[component_stripped] = cbbobj.get_data(text=component, full3D=True)[-1]

        if np.any(water_balance[component_stripped][0] > 0):
            water_balance_summary[component_stripped + '_pos'] = np.sum(
                water_balance[component_stripped][0][water_balance[component_stripped][0] > 0])
        else:
            water_balance_summary[component_stripped + '_pos'] = 0.0
        # End if

        if np.any(water_balance[component_stripped][0] < 0):
            water_balance_summary[component_stripped + '_neg'] = np.sum(
                water_balance[component_stripped][0][water_balance[component_stripped][0] < 0])
        else:
            water_balance_summary[component_stripped + '_neg'] = 0.0
        # End if

        water_bal_summed_titles += component_stripped + '_pos', component_stripped + '_neg'
        water_bal_summed_values += water_balance_summary[component_stripped +
                                                         '_pos'], water_balance_summary[component_stripped + '_neg']
        # End if
    # End for

    water_bal_summed_titles += ['Error']
    Error = sum(water_bal_summed_values)
    water_bal_summed_values += [Error]

    wat_bal_df = pd.DataFrame(water_bal_summed_values, water_bal_summed_titles)
    wat_bal_df.columns = ['Flux m^3/d']
    wat_bal_df = wat_bal_df[wat_bal_df['Flux m^3/d'] != 0.]

    # First step is to set up the plot
    fig = plt.figure(figsize=(20, 10))
    ax = fig.add_subplot(2, 4, 1, aspect='equal')
    ax.set_title('ibound and bc')
    # Next we create an instance of the ModelMap class
    modelmap = flopy.plot.ModelMap(model=self.mf)
    modelmap.plot_ibound()
    bcs = ['RIV', 'GHB']
    self.plot_bcs(modelmap, bcs)

    ax.axes.xaxis.set_ticklabels([])
    ax = fig.add_subplot(2, 4, 2, aspect='equal')
    ax.set_title('Riv flux')
    modelmap = flopy.plot.ModelMap(model=self.mf)
    modelmap.plot_ibound()
    max_val = np.amax(water_balance['RIVER LEAKAGE'].sum(axis=0))
    min_val = np.amin(water_balance['RIVER LEAKAGE'].sum(axis=0))
    if max_val > abs(min_val):
        min_val = -max_val
    else:
        max_val = -min_val
    # End if

    array = modelmap.plot_array(water_balance['RIVER LEAKAGE'].sum(axis=0),
                                alpha=1, cmap='seismic', vmin=min_val, vmax=max_val)
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])
    cbar_ax2 = fig.add_axes([0.43, 0.525, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax2)

    max_val = np.amax(water_balance['HEAD DEP BOUNDS'][1])
    min_val = np.amin(water_balance['HEAD DEP BOUNDS'][1])
    if max_val > abs(min_val):
        min_val = -max_val
    else:
        max_val = -min_val
    # End if

    ax = fig.add_subplot(2, 4, 3, aspect='equal')
    ax.set_title('GHB layer 1')
    modelmap = flopy.plot.ModelMap(model=self.mf)
    modelmap.plot_ibound()
    river_flux = modelmap.plot_array(water_balance['HEAD DEP BOUNDS'][
                                     1], alpha=1, cmap='seismic', vmin=min_val, vmax=max_val)
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])

    cbar_ax1 = fig.add_axes([0.67, 0.525, 0.01, 0.42])
    fig.colorbar(river_flux, cax=cbar_ax1)

    # Setup params to get water balance aspect ratio looking nice
    aspect = float(
        12.5715 / ((wat_bal_df.max()[0] - wat_bal_df.min()[0]) / float(wat_bal_df.shape[1])))

    max_val = np.amax(water_balance['HEAD DEP BOUNDS'][2])
    min_val = np.amin(water_balance['HEAD DEP BOUNDS'][2])
    if max_val > abs(min_val):
        min_val = -max_val
    else:
        max_val = -min_val
    # End if

    ax = fig.add_subplot(2, 4, 4, aspect='equal')
    ax.set_title('GHB layer 2')
    modelmap = flopy.plot.ModelMap(model=self.mf)
    modelmap.plot_ibound()
    recharge = modelmap.plot_array(water_balance['HEAD DEP BOUNDS'][2], masked_values=[
                                   0.], alpha=1, cmap='seismic', vmin=min_val, vmax=max_val)

    start, end = ax.get_xlim()
    start = start // 1000 * 1000 + 1000
    end = end // 1000 * 1000 - 1000
    ax.xaxis.set_ticks(np.arange(start, end, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))

    cbar_ax3 = fig.add_axes([0.19, 0.055, 0.01, 0.42])
    fig.colorbar(recharge, cax=cbar_ax3)

    max_val = np.amax(water_balance['HEAD DEP BOUNDS'][3])
    min_val = np.amin(water_balance['HEAD DEP BOUNDS'][3])
    if max_val > abs(min_val):
        min_val = -max_val
    else:
        max_val = -min_val
    # End if

    ax = fig.add_subplot(2, 4, 5, aspect='equal')
    ax.set_title('GHB layer 3')
    modelmap = flopy.plot.ModelMap(model=self.mf)
    modelmap.plot_ibound()
    elevation = modelmap.plot_array(water_balance['HEAD DEP BOUNDS'][
                                    3], alpha=1, cmap='seismic', vmin=min_val, vmax=max_val)
    ax.yaxis.set_ticklabels([])
    ax.xaxis.set_ticks(np.arange(start, end, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))

    max_val = np.amax(water_balance['HEAD DEP BOUNDS'][4])
    min_val = np.amin(water_balance['HEAD DEP BOUNDS'][4])
    if max_val > abs(min_val):
        min_val = -max_val
    else:
        max_val = -min_val
    # End if

    ax = fig.add_subplot(2, 4, 6, aspect='equal')
    ax.set_title('GHB layer 4')
    modelmap = flopy.plot.ModelMap(model=self.mf)
    modelmap.plot_ibound()
    elevation = modelmap.plot_array(water_balance['HEAD DEP BOUNDS'][
                                    4], alpha=1, cmap='seismic', vmin=min_val, vmax=max_val)
    ax.yaxis.set_ticklabels([])
    ax.xaxis.set_ticks(np.arange(start, end, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    
    
    cbar_ax4 = fig.add_axes([0.43, 0.055, 0.01, 0.42])
    fig.colorbar(elevation, cax=cbar_ax4)

    max_val = np.amax(water_balance['HEAD DEP BOUNDS'][5])
    min_val = np.amin(water_balance['HEAD DEP BOUNDS'][5])
    if max_val > abs(min_val):
        min_val = -max_val
    else:
        max_val = -min_val
    # End if

    ax = fig.add_subplot(2, 4, 7, aspect='equal')
    ax.set_title('GHB layer 5')
    modelmap = flopy.plot.ModelMap(model=self.mf)
    modelmap.plot_ibound()
    storage = modelmap.plot_array(water_balance['HEAD DEP BOUNDS'][
                                  5], alpha=1, cmap='seismic', vmin=min_val, vmax=max_val)
    ax.yaxis.set_ticklabels([])
    ax.xaxis.set_ticks(np.arange(start, end, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))

    cbar_ax5 = fig.add_axes([0.67, 0.055, 0.01, 0.42])
    fig.colorbar(storage, cax=cbar_ax5)

    max_val = np.amax(water_balance['HEAD DEP BOUNDS'][6])
    min_val = np.amin(water_balance['HEAD DEP BOUNDS'][6])
    if max_val > abs(min_val):
        min_val = -max_val
    else:
        max_val = -min_val
    # End if

    ax = fig.add_subplot(2, 4, 8, aspect='equal')
    ax.set_title('GHB layer 6')
    modelmap = flopy.plot.ModelMap(model=self.mf)
    modelmap.plot_ibound()
    storage = modelmap.plot_array(water_balance['HEAD DEP BOUNDS'][
                                  6], alpha=1, cmap='seismic', vmin=min_val, vmax=max_val)
    ax.yaxis.set_ticklabels([])
    ax.xaxis.set_ticks(np.arange(start, end, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))

    cbar_ax5 = fig.add_axes([0.67, 0.055, 0.01, 0.42])
    fig.colorbar(storage, cax=cbar_ax5)
    
    
    fig.subplots_adjust(left=0.01, right=0.95, bottom=0.05, top=0.95,
                        wspace=0.1, hspace=0.12)
    plt.show()

    cbbobj.close()
    headobj.close()
# End viewGHB()

def plot_bc_by_zone(self, bc_name, nper=-1):
    """TODO: Docs"""

    cbbobj = bf.CellBudgetFile(os.path.join(self.data_folder, self.name + '.cbc'))

    water_balance_components = cbbobj.textlist
    water_balance = {}
    water_balance_summary = {}
    water_bal_summed_titles = []
    water_bal_summed_values = []

    for component in water_balance_components:
        component_stripped = component.lstrip().rstrip()
        if 'FLOW' in component_stripped:
            continue
        # End if
        water_balance[component_stripped] = cbbobj.get_data(text=component, full3D=True)[nper]

        if np.any(water_balance[component_stripped][0] > 0):
            water_balance_summary[component_stripped + '_pos'] = np.sum(
                water_balance[component_stripped][0][water_balance[component_stripped][0] > 0])
        else:
            water_balance_summary[component_stripped + '_pos'] = 0.0
        # End if

        if np.any(water_balance[component_stripped][0] < 0):
            water_balance_summary[component_stripped + '_neg'] = np.sum(
                water_balance[component_stripped][0][water_balance[component_stripped][0] < 0])
        else:
            water_balance_summary[component_stripped + '_neg'] = 0.0
        # End if

        water_bal_summed_titles += component_stripped + '_pos', component_stripped + '_neg'
        water_bal_summed_values += water_balance_summary[component_stripped +
                                                         '_pos'], water_balance_summary[component_stripped + '_neg']
        # End if
    # End for

    water_bal_summed_titles += ['Error']
    Error = sum(water_bal_summed_values)
    water_bal_summed_values += [Error]

    wat_bal_df = pd.DataFrame(water_bal_summed_values, water_bal_summed_titles)
    wat_bal_df.columns = ['Flux m^3/d']
    wat_bal_df = wat_bal_df[wat_bal_df['Flux m^3/d'] != 0.]

    wb_bc = water_balance[bc_name]
    
    if nper == 'all':
        wb_bc = np.mean(wb_bc, axis=0)
        zoned = self.output_var_by_zone(wb_bc)
        wb_bc = zoned
    else:
        zoned = self.output_var_by_zone(wb_bc)
        wb_bc = zoned
    # End if

    cbbobj.close()

    # First step is to set up the plot
    fig = plt.figure(figsize=(20, 10))

    layers = {0: 'Coonambidgal',
              1: 'Basalt',
              2: 'Shepparton',
              3: 'utam',
              4: 'Calivil',
              5: 'Renmark',
              6: 'Basement'}

    for key in layers:          
        ax = fig.add_subplot(2, 4, key + 1, aspect='equal')
        ax.set_title(layers[key])
        # Next we create an instance of the ModelMap class
        modelmap = flopy.plot.ModelMap(model=self.mf)
        modelmap.plot_ibound()

        array = modelmap.plot_array(wb_bc[key], masked_values= [0.],
                                    alpha=1, cmap='seismic')
        
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        cbar_ax2 = fig.add_axes([0.95, 0.98, 0.1, 0.42])
        fig.colorbar(array, cax=cbar_ax2)
    
 # End plot_bc_by_zone()
