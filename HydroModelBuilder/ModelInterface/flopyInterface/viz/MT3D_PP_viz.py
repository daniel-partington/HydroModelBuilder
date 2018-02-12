import os

import flopy
import matplotlib.pyplot as plt
import numpy as np


def compareAllObs(self):
    """TODO: Docs"""

    concobj = self.import_concs()
    times = concobj.get_times()

    scatterx = []
    scattery = []
    obs_sim_zone_all = []

    # The definition of obs_sim_zone looks like:
    for i in range(self.mf_model.model_data.model_time.t['steps']):
        conc = concobj.get_data(totim=times[i])
        self.compare_observed('C14', conc, nper=i)
        obs_sim_zone_all += self.obs_sim_zone

    scatterx = np.array([h[0] for h in obs_sim_zone_all])
    scattery = np.array([h[1] for h in obs_sim_zone_all])

    # First step is to set up the plot
    width = 20
    height = 5
    multiplier = 1.0
    fig = plt.figure(figsize=(width * multiplier, height * multiplier))

    ax = fig.add_subplot(1, 3, 1)
    ax.set_title('Residuals')
    ax.hist([loc[0] - loc[1] for loc in obs_sim_zone_all], bins=20, alpha=0.5)

    ax = fig.add_subplot(1, 3, 2)
    ax.set_title('Sim vs Obs (%d points)' % (len(scatterx)))

    comp_zone_plots = {}
    colours = ['r', 'orangered', 'y', 'green', 'teal', 'blue', 'fuchsia']
    for i in xrange(1, 8):
        scatterx2 = [loc[0] for loc in obs_sim_zone_all if loc[2] == float(i)]
        scattery2 = [loc[1] for loc in obs_sim_zone_all if loc[2] == float(i)]
        # print len(scatterx2), colours[i-1]
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
            sum1 += num1 ** np.float64(2.)
            sum2 += num2 ** np.float64(2.)

        ME = 1 - sum1 / sum2

        ax.text(150, 75, 'Model Efficiency = %4.2f' % (ME))

        # for PBIAS
        def pbias(simulated, observed):
            """
            :param simulated:
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

    ax = fig.add_subplot(1, 3, 3)
    ax.set_title('Residuals in space')

    modelmap = flopy.plot.ModelMap(model=self.mf_model.mf)
    modelmap.plot_ibound()

    x = np.array([h[3] for h in obs_sim_zone_all])
    y = np.array([h[4] for h in obs_sim_zone_all])
    zone = [h[2] for h in obs_sim_zone_all]
    residuals = [h[0] - h[1] for h in obs_sim_zone_all]
    residuals = np.absolute(residuals)

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
                rgb_ref += [rgb_all[index]]
            # End if
        # End for
    # End for

    zone = np.array(zone)
    rgba_colors = np.zeros((len(x), 4))
    # for red the first column needs to be one
    for i in range(1, 8):
        rgba_colors[:, 0][zone == i] = rgb_ref[i - 1][0]
        rgba_colors[:, 1][zone == i] = rgb_ref[i - 1][1]
        rgba_colors[:, 2][zone == i] = rgb_ref[i - 1][2]
    # the fourth column needs to be your alphas
    rgba_colors[:, 3] = residuals / np.max(residuals)  # alphas

    plt.scatter(x, y, color=rgba_colors)
    plt.show()
# End compareAllObs()


def viewConcsByZone(self, nper='all', specimen=None):
    """
    :param nper: (Default value = 'all')

    :param specimen: (Default value = None)
    """

    # Create the headfile object
    concobj = self.import_concs()
    times = concobj.get_times()
    if nper == 'all':
        conc = concobj.get_alldata()
        conc = np.mean(conc, axis=0)
        zoned = self.ConcsByZone(conc)
        conc = zoned
    elif nper == 'final':
        conc = concobj.get_data(totim=times[-1])
        zoned = self.ConcsByZone(conc)
        conc = zoned
    else:
        conc = concobj.get_data(totim=times[nper])
        zoned = self.ConcsByZone(conc)
        conc = zoned
    # End if

    # First step is to set up the plot
    width = 20
    height = 10
    multiplier = 1.
    fig = plt.figure(figsize=(width * multiplier, height * multiplier))

    vmin = np.amin(conc[conc > 0.])
    vmax = np.amax(conc)

    ax = fig.add_subplot(2, 4, 1, aspect='equal')

    ax.set_title('ibound and bc')
    # Next we create an instance of the ModelMap class
    modelmap = flopy.plot.ModelMap(model=self.mf_model.mf)
    modelmap.plot_ibound()

    modelmap.plot_bc('RIV', plotAll=True)
    try:
        modelmap.plot_bc('WEL', plotAll=True)
    except Exception:
        pass
    modelmap.plot_bc('GHB', plotAll=True)
    modelmap.plot_bc('SFR', plotAll=True)
    try:
        modelmap.plot_bc('DRN', plotAll=True)
    except Exception:
        pass
    ax.axes.xaxis.set_ticklabels([])

    ax = fig.add_subplot(2, 4, 2, aspect='equal')
    ax.set_title('Coonambidgal')
    modelmap = flopy.plot.ModelMap(model=self.mf_model.mf)
    min_conc = -100.0
    max_conc = 100.0
    temp = max_conc
    max_conc = vmax
    vmax = 100.0

    array = modelmap.plot_array(
        conc[0], masked_values=[-999.98999023, max_conc, min_conc], alpha=0.5, vmin=vmin, vmax=vmax)

    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])
    cbar_ax2 = fig.add_axes([0.43, 0.525, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax2)

    ax = fig.add_subplot(2, 4, 3, aspect='equal')
    ax.set_title('Shepparton')
    modelmap = flopy.plot.ModelMap(model=self.mf_model.mf)
    array = modelmap.plot_array(
        conc[2], masked_values=[-999.98999023, max_conc, min_conc, np.nan], alpha=0.5, vmin=vmin, vmax=vmax)
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])
    cbar_ax1 = fig.add_axes([0.67, 0.525, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax1)

    ax = fig.add_subplot(2, 4, 5, aspect='equal')
    ax.set_title('Calivil')
    # , sr=self.mf.dis.sr, dis=self.mf.dis)
    modelmap = flopy.plot.ModelMap(model=self.mf_model.mf)
    array = modelmap.plot_array(
        conc[4], masked_values=[-999.98999023, max_conc, min_conc, np.nan], alpha=0.5, vmin=vmin, vmax=vmax)
    start, end = ax.get_xlim()
    start = start // 1000 * 1000 + 1000
    end = end // 1000 * 1000 - 1000
    ax.xaxis.set_ticks(np.arange(start, end, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    cbar_ax3 = fig.add_axes([0.19, 0.055, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax3)

    ax = fig.add_subplot(2, 4, 6, aspect='equal')
    ax.set_title('Renmark')
    # , sr=self.mf.dis.sr, dis=self.mf.dis)
    modelmap = flopy.plot.ModelMap(model=self.mf_model.mf)
    array = modelmap.plot_array(
        conc[5], masked_values=[-999.98999023, max_conc, min_conc, np.nan], alpha=0.5, vmin=vmin, vmax=vmax)
    ax.yaxis.set_ticklabels([])
    ax.xaxis.set_ticks(np.arange(start, end, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    cbar_ax4 = fig.add_axes([0.43, 0.055, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax4)

    ax = fig.add_subplot(2, 4, 7, aspect='equal')
    ax.set_title('Basement')
    modelmap = flopy.plot.ModelMap(model=self.mf_model.mf)
    array = modelmap.plot_array(
        conc[6], masked_values=[-999.98999023, max_conc, min_conc, np.nan], alpha=0.5, vmin=vmin, vmax=vmax)
    ax.yaxis.set_ticklabels([])
    ax.xaxis.set_ticks(np.arange(start, end, 20000.0))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    cbar_ax5 = fig.add_axes([0.67, 0.055, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax5)
    fig.subplots_adjust(left=0.01, right=0.95, bottom=0.05, top=0.95, wspace=0.1, hspace=0.12)

    plt.show()
# End viewConcsByZone()


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

    #colours = ['r', 'orangered', 'y', 'green', 'teal', 'blue', 'fuchsia']
    ax.hist(residuals, bins=20, alpha=0.5, color='black', histtype='step', label='all')

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
            return np.sum(simulated - observed) * 100 / np.sum(observed)

        ax.text(xmin + 0.45 * (xmax - xmin), ymin + 0.3 * (ymax - ymin),
                'PBIAS = %4.2f%%' % (pbias(scattery, scatterx)))

        # For rmse
        def rmse(simulated, observed):
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


def compareAllObs2(self, specimen):

    conc = None
    sft_conc = None
    obs_group = self.mf_model.model_data.observations.obs_group

    obs_sim_zone_all = []

    # Write observation to file
    for obs_set in obs_group:

        obs_sim_zone_all = []

        obs_type = obs_group[obs_set]['obs_type']
        # Import the required model outputs for processing
        if obs_type not in ['concentration', 'EC', 'Radon']:
            continue
        else:
            print("Processing {}".format(obs_set))
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

        if obs_group[obs_set]['domain'] == 'stream':
            sft_location = obs_group[obs_set]['locations']['seg_loc']

            for observation in obs_df.index:
                interval = int(obs_df['interval'].loc[observation])
                name = obs_df['name'].loc[observation]
                obs = obs_df['value'].loc[observation]
                seg = sft_location.loc[name]
                sft = sft_conc
                times = sft['TIME'].unique()
                col_of_interest = obs_type
                if obs_type == 'EC':
                    col_of_interest = 'SFR-CONCENTRATION'
                sim_obs = sft[(sft['SFR-NODE'] == seg) &
                              (sft['TIME'] == times[interval])][col_of_interest].tolist()[0]
                obs_sim_zone_all += [[obs, sim_obs, seg]]

        if obs_group[obs_set]['domain'] == 'porous':
            for observation in obs_df.index:
                interval = int(obs_df['interval'].loc[observation])
                name = obs_df['name'].loc[observation]
                obs = obs_df['value'].loc[observation]
                zone = obs_df['zone'].loc[observation]

                (lay, row, col) = [sim_map_dict[name][0],
                                   sim_map_dict[name][1], sim_map_dict[name][2]]
                sim_conc = [conc[interval][lay][row][col]]
                sim_obs = np.mean(sim_conc)
                obs_sim_zone_all += [[obs, sim_obs, zone]]
            # End for
        # End if

        self._plot_obs_vs_sim(obs_set, obs_sim_zone_all, unc=2)


# End compareAllObs()
