import os

import flopy
import flopy.utils.binaryfile as bf
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import six
from flopy.utils.sfroutputfile import SfrFile
from matplotlib import colors


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


def SFRoutput_plot(self):
    """TODO: Docs"""

    raise NotImplementedError("Tried to call unfinished method!")
    sfrout = SfrFile(os.path.join(self.data_folder, self.name + ".sfr.out"))
    sfr_df = sfrout.get_dataframe()
    sfr_df
# End SFRoutput_plot()


def compareAllObs_metrics(self, to_file=False):

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

    scatterx = np.array([h[0] for h in obs_sim_zone_all])
    scattery = np.array([h[1] for h in obs_sim_zone_all])

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

        # for PBIAS
        def pbias(simulated, observed):
            return np.sum(simulated - observed) * 100 / np.sum(observed)

        PBIAS = pbias(scattery, scatterx)

        # For rmse
        def rmse(simulated, observed):
            return np.sqrt(((simulated - observed) ** 2).mean())

        RMSE = rmse(scattery, scatterx)

    if to_file:
        with open(os.path.join(self.data_folder, 'Head_Obs_Model_Measures.txt'), 'w') as f:
            f.write('ME PBIAS RMSE\n')
            f.write('{} {} {}'.format(ME, PBIAS, RMSE))

    return ME, PBIAS, RMSE
# End compareAllObs_metrics()


def compareAllObs(self):
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
        self.compare_observed('head', head, nper=i)
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

    sum1 = 0.0
    sum2 = 0.0

    if len(scatterx) != 0:
        mean = np.mean(scatterx)
        for i in range(len(scatterx)):
            num1 = (scatterx[i] - scattery[i])
            num2 = (scatterx[i] - mean)
            sum1 += num1**np.float64(2.0)
            sum2 += num2**np.float64(2.0)

        ME = 1 - sum1 / sum2

        ymin, ymax = ax.get_ylim()
        xmin, xmax = ax.get_xlim()

        ax.text(xmin + 0.45 * (xmax - xmin), ymin + 0.4 * (ymax - ymin),
                'Model Efficiency = %4.2f' % (ME))

        # for PBIAS
        def pbias(simulated, observed):
            """
            :param simulated: param observed:

            :param observed:
            """

            return np.sum(simulated - observed) * 100 / np.sum(observed)

        ax.text(xmin + 0.45 * (xmax - xmin), ymin + 0.3 * (ymax - ymin),
                'PBIAS = %4.2f%%' % (pbias(scattery, scatterx)))

        # For rmse
        def rmse(simulated, observed):
            """
            :param simulated: param observed:

            :param observed:
            """

            return np.sqrt(((simulated - observed) ** 2).mean())

        ax.text(xmin + 0.45 * (xmax - xmin), ymin + 0.2 * (ymax - ymin),
                'RMSE = %4.2f' % (rmse(scattery, scatterx)))

    ax.plot(ax.get_ylim(), ax.get_ylim())

    ax = fig.add_subplot(1, 3, 3)
    ax.set_title('Residuals in space')

    modelmap = flopy.plot.ModelMap(model=self.mf)
    modelmap.plot_ibound()

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

    plt.colorbar()
    plt.show()
# End compareAllObs()


def _plot_obs_vs_sim(self, obs_set, obs_sim_zone_all, unc=None):
    scatterx = np.array([h[0] for h in obs_sim_zone_all])
    scattery = np.array([h[1] for h in obs_sim_zone_all])

    residuals = [loc[0] - loc[1] for loc in obs_sim_zone_all]

    # First step is to set up the plot
    width = 20
    height = 5
    multiplier = 1.
    fig = plt.figure(figsize=(width * multiplier, height * multiplier))

    ax = fig.add_subplot(1, 2, 1)  # , aspect='equal')
    ax.set_title('Residuals')

    ax.hist(residuals, bins=20, alpha=0.5, color='black', histtype='step', label='all')
    
    for i in range(self.model_data.model_time.t['steps']):
        head = headobj.get_data(totim=times[i])
        self.compare_observed('head', head, nper=i)
        obs_sim_zone_all += self.obs_sim_zone

    plt.legend(loc='upper left', ncol=4, fontsize=11)

    ax = fig.add_subplot(1, 2, 2)
    ax.set_title('{}: Sim vs Obs ({} points)'.format(obs_set.upper(), len(scatterx)))
    ax.scatter(scatterx, scattery, facecolors='none', alpha=0.5)

    plt.xlabel('Observed')
    plt.ylabel('Simulated', labelpad=10)

    sum1 = 0.0
    sum2 = 0.0

    if len(scatterx) != 0:
        mean = np.mean(scatterx)
        for i in range(len(scatterx)):
            num1 = (scatterx[i] - scattery[i])
            num2 = (scatterx[i] - mean)
            sum1 += num1**np.float64(2.0)
            sum2 += num2**np.float64(2.0)

        ME = 1 - sum1 / sum2

        ymin, ymax = ax.get_ylim()
        xmin, xmax = ax.get_xlim()

        ax.text(xmin + 0.45 * (xmax - xmin), ymin + 0.4 * (ymax - ymin),
                'Model Efficiency = %4.2f' % (ME))

        # for PBIAS
        def pbias(simulated, observed):
            """
            :param simulated: param observed:

            :param observed:
            """

            return np.sum(simulated - observed) * 100 / np.sum(observed)

        ax.text(xmin + 0.45 * (xmax - xmin), ymin + 0.3 * (ymax - ymin),
                'PBIAS = %4.2f%%' % (pbias(scattery, scatterx)))

        # For rmse
        def rmse(simulated, observed):
            """
            :param simulated: param observed:

            :param observed:
            """

            return np.sqrt(((simulated - observed) ** 2).mean())

        ax.text(xmin + 0.45 * (xmax - xmin), ymin + 0.2 * (ymax - ymin),
                'RMSE = %4.2f' % (rmse(scattery, scatterx)))

    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    new = (min(xlim[0], ylim[0]), max(xlim[1], ylim[1]))
    new_upper = (min(xlim[0], ylim[0]) + unc, max(xlim[1], ylim[1]) + unc)
    new_lower = (min(xlim[0], ylim[0]) - unc, max(xlim[1], ylim[1]) - unc)
    ax.plot(new, new, color='grey')
    ax.plot(new_upper, new_upper, color='grey')
    ax.plot(new_lower, new_lower, color='grey')
    ax.fill_between(new, new_lower, new_upper, color='grey', alpha=0.3)
    ax.set_xlim(new)
    ax.set_ylim(new)


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
                    zone_txt = 'head'
                else:
                    zone_txt = zone
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

        self._plot_obs_vs_sim(obs_set, obs_sim_zone_all, unc=2)


# End compareAllObs()

def viewHeadsByZone(self, nper='all'):
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
        self.compare_observed('head', head_orig, nper=nper)
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

    modelmap.plot_bc('RIV')  # Can also plot 'WEL' and 'DRN'
    modelmap.plot_bc('SFR')
    modelmap.plot_bc('GHB')
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


def viewHeadsByZone2(self, iter_num, nper='all'):
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
            self.compare_observed('head', head_orig, nper=i)
            scatterx += [h[0] for h in self.obs_sim_zone]
            scattery += [h[1] for h in self.obs_sim_zone]
            obs_sim_zone_all += self.obs_sim_zone
        self.obs_sim_zone = obs_sim_zone_all
    else:
        self.compare_observed('head', head_orig, nper=nper)
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

    modelmap.plot_bc('RIV', plotAll=True)
    modelmap.plot_bc('WEL', plotAll=True)
    modelmap.plot_bc('GHB', plotAll=True)
    try:
        modelmap.plot_bc('SFR', plotAll=True)
    except Exception:
        print("No SFR package present")
    # End try

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

        ax.text(150, 75, 'Model Efficiency = %4.2f' % (ME))

        # for PBIAS
        def pbias(simulated, observed):
            """
            :param simulated: param observed:

            :param observed:
            """

            return np.sum(simulated - observed) * 100 / np.sum(observed)

        ax.text(150, 40, 'PBIAS = %4.2f%%' % (pbias(scattery, scatterx)))

        # For rmse
        def rmse(simulated, observed):
            """
            :param simulated: param observed:

            :param observed:
            """

            return np.sqrt(((simulated - observed) ** 2).mean())

        ax.text(150, 20, 'RMSE = %4.2f' % (rmse(scattery, scatterx)))

    ax.plot(ax.get_ylim(), ax.get_ylim())
    fig.subplots_adjust(left=0.01, right=0.95, bottom=0.05, top=0.95, wspace=0.1, hspace=0.12)
    plt.savefig('run_viewheads_{}.png'.format(iter_num), bbox_inches='tight')

# End viewHeadsByZone2()


def viewHeads(self):
    """TODO: Docs"""

    # Create the headfile object
    headobj = self.importHeads()
    cbbobj = self.importCbb()
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
    modelmap.plot_bc('RIV', plotAll=True)
    modelmap.plot_bc('SFR', plotAll=True)
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
    modelmap.plot_bc('RIV')
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

    modelmap.plot_bc('RIV', alpha=0.3)
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
    modelmap.plot_bc('RIV', plotAll=True)
    modelmap.plot_bc('SFR', plotAll=True)
    try:
        modelmap.plot_bc('GHB', plotAll=True)
    except Exception:
        pass
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
    river_flux = modelmap.plot_array(water_balance['STREAM LEAKAGE'].sum(axis=0), alpha=0.5)
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
    elev = np.ma.masked_where(self.model_data.model_mesh3D[1][0] == 0, elev)
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
        pressure, masked_values=[x for x in np.unique(pressure) if x > 0.], alpha=0.5)  # , vmin=-50, vmax=50)
    self.plotRiverFromRiverSegData(ax)
    ax.yaxis.set_ticklabels([])
    ax.xaxis.set_ticks(np.arange(start, end, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))

    cbar_ax5 = fig.add_axes([0.67, 0.055, 0.01, 0.42])
    fig.colorbar(storage, cax=cbar_ax5)

    fig.subplots_adjust(left=0.01, right=0.95, bottom=0.05, top=0.95,
                        wspace=0.1, hspace=0.12)
    plt.show()

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
    modelmap.plot_bc('RIV')

    try:
        modelmap.plot_bc('GHB')
    except Exception as e:
        print(e)
    # End try

    ax.axes.xaxis.set_ticklabels([])
    ax = fig.add_subplot(2, 4, 2, aspect='equal')
    ax.set_title('Riv flux')
    modelmap = flopy.plot.ModelMap(model=self.mf)
    modelmap.plot_ibound()
    max_val = np.amax(water_balance['RIVER LEAKAGE'][0])
    min_val = np.amin(water_balance['RIVER LEAKAGE'][0])
    if max_val > abs(min_val):
        min_val = -max_val
    else:
        max_val = -min_val
    # End if

    array = modelmap.plot_array(water_balance['RIVER LEAKAGE'][
                                0], alpha=1, cmap='seismic', vmin=min_val, vmax=max_val)
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
    ax.set_title('GHB layer 2')
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

    ax = fig.add_subplot(2, 4, 5, aspect='equal')
    ax.set_title('GHB layer 3')
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

    ax = fig.add_subplot(2, 4, 6, aspect='equal')
    ax.set_title('GHB layer 4')
    modelmap = flopy.plot.ModelMap(model=self.mf)
    modelmap.plot_ibound()
    elevation = modelmap.plot_array(water_balance['HEAD DEP BOUNDS'][
                                    3], alpha=1, cmap='seismic', vmin=min_val, vmax=max_val)
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

    fig.subplots_adjust(left=0.01, right=0.95, bottom=0.05, top=0.95,
                        wspace=0.1, hspace=0.12)
    plt.show()
# End viewGHB()
