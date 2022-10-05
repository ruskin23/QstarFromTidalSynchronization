from fileinput import filename
from turtle import color
import emcee
import sys
from configargparse import ArgumentParser
import corner
import matplotlib.pyplot as plt
import numpy
import h5py
import glob

from pathlib import Path
from directories import directories


home_dir=str(Path.home())
path=directories(home_dir)
sys.path.append(path.poet_path+'/PythonPackage')
sys.path.append(path.poet_path+'/scripts')
sys.path.append('/home/ruskin/projects/')
sys.path.append('/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/version2_emcee/h5analysis/')
sys.path.append('/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/version2_emcee/')


from orbital_evolution.transformations import phase_lag, lgQ
from general_purpose_python_modules.emcee_quantile_convergence import *
from utils import *
from h5_analysis import *
import scipy

_quantities=[scipy.stats.norm.cdf(c) for c in [-2,-1,1,2]]

def get_burnin(system_kic):

    burnins = []
    with open('period_dependence/{}.txt'.format(system_kic),'r') as f:
        next(f)
        for lines in f:
            x = lines.split()
            for i in [3,6,9,12]:
                burnins.append(int(x[i].split('/')[0]))
                max_step = int(x[i].split('/')[1])

    max_burnins = max(burnins)
    return max_burnins


def cut_off_periods(system_kic):

    differences = []
    periods = []
    
    data = numpy.empty((50,5),float)
    with open('period_dependence/{}.txt'.format(system_kic),'r') as f:
        next(f)
        for i,lines in enumerate(f):
            x = lines.split()
            differences.append(abs(float(x[1]) - float(x[10])))
            periods.append(float(x[0]))

            data[i] = numpy.array([[float(x[k]) for k in [0,1,4,7,10]]])

    min_index = differences.index(min(differences))
    #print(min_index)
    if min_index < 3:
        period_cut_off = periods[:min_index+3]
    elif min_index > len(periods) - 4:
        period_cut_off = periods[min_index-3:]
    else:
        period_cut_off = periods[min_index-3:min_index+3]
    
    #print(system_kic, period_cut_off)
    # data = numpy.transpose(data)
    # for k in range(4):
    #     plt.plot(data[0,:i],data[k+1,:i])
    # plt.vlines(x = period_cut_off[0], color = 'b', ymin = 5, ymax = 12)
    # plt.vlines(x = period_cut_off[-1], color = 'b', ymin = 5, ymax = 12)
    # plt.show()
    
    return (period_cut_off[0],period_cut_off[1])


def plot_distribution(system_filename):

    system_kic = system_filename.split('.')[0].split('_')[1]

    prior_samples = corrected_samples(system_filename)
    period_min,period_max = cut_off_periods(system_kic)

    lgQ_array = numpy.linspace(5,20,500)
    period_array = numpy.linspace(numpy.log10(0.5), numpy.log10(50), 500)
    data = numpy.empty((500, 500), float)

    for i,periods in enumerate(period_array):
        if periods < period_min or periods > period_max:
            data[i] = numpy.ones(500)/500.0
        else:
            #print(i)
            lgQ_samples = get_dissipation_samples(prior_samples, tidal_period=periods)
            burn_in = get_burnin(system_kic)
            #print(burn_in)
            #print(len(lgQ_samples))
            #print(lgQ_samples)
            lgQ_samples = lgQ_samples[burn_in+1:]
            #print(lgQ_samples)
            lgQ_samples = lgQ_samples.flatten()
            #print(lgQ_samples)
            _dist = DiscreteSampling(lgQ_samples, 0.1)
            # for q in _quantities:
            #     plt.scatter(periods,_dist._ppf(q))
            #     print(_dist._ppf(q))

            data[i] = _dist._pdf(lgQ_array)
            data[i] = data[i]/numpy.sum(data[i])
    
    # converged_data = numpy.empty((50,5), float)
    # with open(f'period_dependence/{system_kic}.txt','r') as f:
    #     next(f)
    #     for i,lines in enumerate(f):
    #         x = lines.split()
    #         converged_data[i] = numpy.array([float(x[k]) for k in [0,1,4,7,10]])

    # converged_data = numpy.transpose(converged_data)
    # for k in range(4):
    #     plt.plot(converged_data[0,:i],converged_data[k+1,:i])
    # plt.show()

    return data


def joint_distribution():

    with open('period_dependence/stop_systems.txt','r') as f:
        system_kic = f.read().split('\n')
        system_kic = [sk for sk in system_kic if sk != '']
    
    # system_kic = ['7987749']
    #print(system_kic)
    data = numpy.ones((500, 500))
    for S in system_kic:
        #print('\nSystem  = {}'.format(S))
        filename = f'system_{S}.h5'

        d_i = plot_distribution(filename)
        data = data*d_i

    period_array = numpy.linspace(numpy.log10(0.5), numpy.log10(50), 500)
    lgQ_array = numpy.linspace(5,12,500)
    
    quantile_grid = numpy.empty((500,4), float)
    ignore_samples = []
    for i, dist in enumerate(data):
        normalized_dist = dist/numpy.sum(dist)
        if normalized_dist[0] == normalized_dist[-1]:
            ignore_samples.append(period_array[i])
        #print(normalized_dist)
        combined_dist = DiscreteDistributions(normalized_dist, lgQ_array)
        quantile_grid[i] = numpy.array([combined_dist._ppf(q) for q in _quantities])

    # quantile_grid = numpy.transpose(quantile_grid)
    # for i in range(4):
    #     plt.plot(period_array, quantile_grid[i])
    minQ2_array = []
    maxQ2_array = []
    minQ1_array = []
    maxQ1_array = []
    p_reduced = []
    for q, p in zip(quantile_grid, period_array):
        minQ2 = q[0]
        maxQ2 = q[3]
        minQ1 = q[1]
        maxQ1 = q[2]
        if p not in ignore_samples:
            #print(p)
            p_reduced.append(p)
            minQ2_array.append(minQ2)
            maxQ2_array.append(maxQ2)
            minQ1_array.append(minQ1)
            maxQ1_array.append(maxQ1)

            plt.vlines(p, ymin=minQ2, ymax=maxQ2, linewidth=0.5)
            plt.vlines(p, ymin=minQ1, ymax=maxQ1, linewidth=1.2)
    plt.savefig('a1.png')
    plt.close()

    # quantile_grid = numpy.transpose(quantile_grid)
    minQ2 = q[0]
    maxQ2 = q[3]
    minQ1 = q[1]
    maxQ1 = q[2]
    plt.fill_between(p_reduced, minQ2_array, minQ1_array, color='green', alpha=0.2)
    plt.fill_between(p_reduced, minQ1_array, maxQ1_array, color='red', alpha=0.5)
    plt.fill_between(p_reduced, maxQ1_array, maxQ2_array, color='green', alpha=0.2)
    plt.savefig('a2.png')



if __name__ == '__main__':

    args = cmdline_args()
    system_filename = args.chains_file
    
    # plot_distribution(system_filename)
    joint_distribution()
