import numpy as np
import matplotlib.pyplot as plt

# Model Efficiency
def metric_me(simulated, observed):
    '''
    Calculates the model efficiency. Length of both argument lists should be
    equal
    
    :param simulated: list[float], list of simulated observations
    :param observed: list[float], list of observed observations
    
    :returns: float    
    '''
    sum1 = 0.0
    sum2 = 0.0

    mean = np.mean(observed)
    for i in range(len(observed)):
        num1 = (observed[i] - simulated[i])
        num2 = (simulated[i] - mean)
        sum1 += num1 ** np.float64(2.)
        sum2 += num2 ** np.float64(2.)

    return 1 - sum1 / sum2

# Percent BIAS
def metric_pbias(simulated, observed):
    '''
    Calculates the percent bias. Length of both argument lists should be
    equal
    
    :param simulated: list[float], list of simulated observations
    :param observed: list[float], list of observed observations
    
    :returns: float    
    '''
    return np.sum(simulated - observed) * 100 / np.sum(observed)

# For root mean squared error
def metric_rmse(simulated, observed):
    '''
    Calculates the root mean square error. Length of both argument lists should be
    equal
    
    :param simulated: list[float], list of simulated observations
    :param observed: list[float], list of observed observations
    
    :returns: float    
    '''
    return np.sqrt(((simulated - observed) ** 2).mean())
    
def plot_obs_vs_sim(obs_set, obs_sim_zone_all, unc=None):
    '''Plot of observed vs simulated
    
    :param obs_set: str, name of the observation set
    :param obs_sim_zone_all: list[list[float,float]], list of lists of observed and simulated
    :param unc: float, uncertainty associated with the observation type which is used to plot
                grey area around the centre line of perfect fit
    
    :returns: matplotlib axis object 
    '''
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

    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()

    ax.text(xmin + 0.45 * (xmax - xmin), ymin + 0.4 * (ymax - ymin),
            'Model Efficiency = %4.2f' % (metric_me(scattery, scatterx)))

    ax.text(xmin + 0.45 * (xmax - xmin), ymin + 0.3 * (ymax - ymin),
            'PBIAS = %4.2f%%' % (metric_pbias(scattery, scatterx)))

    ax.text(xmin + 0.45 * (xmax - xmin), ymin + 0.2 * (ymax - ymin),
            'RMSE = %4.2f' % (metric_rmse(scattery, scatterx)))

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
    
    return ax

