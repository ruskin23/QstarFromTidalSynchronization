import pickle
import dill
import numpy as np

import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from matplotlib import rcParams
rcParams.update({'xtick.major.pad': '7.0'})
rcParams.update({'xtick.major.size': '7.5'})
rcParams.update({'xtick.major.width': '1.5'})
rcParams.update({'xtick.minor.pad': '7.0'})
rcParams.update({'xtick.minor.size': '3.5'})
rcParams.update({'xtick.minor.width': '1.0'})
rcParams.update({'ytick.major.pad': '7.0'})
rcParams.update({'ytick.major.size': '7.5'})
rcParams.update({'ytick.major.width': '1.5'})
rcParams.update({'ytick.minor.pad': '7.0'})
rcParams.update({'ytick.minor.size': '3.5'})
rcParams.update({'ytick.minor.width': '1.0'})
rcParams.update({'font.size': 30})

import dynesty
from dynesty import plotting as dyplot

def result_plots():
    # 3-D plots of position and likelihood, colored by weight
    fig = plt.figure(figsize=(30, 10))
    ax = fig.add_subplot(121, projection='3d')

    # plotting the initial run
    p = ax.scatter(custom_results.samples[:, 0], custom_results.samples[:, 1],  custom_results.samples[:, 2],
                   marker='o', c=np.exp(custom_results.logwt) * 1e7, linewidths=(0.,), cmap='coolwarm')
    #ax.set_xlim(-10., 10.)
    #ax.set_xticks(np.linspace(-10., 10., 5))
    #ax.set_xlabel(r'$x$', labelpad=25)
    #ax.set_ylim(-10., 10.)
    #ax.set_yticks(np.linspace(-10., 10., 5))
    #ax.set_ylabel(r'$y$', labelpad=25)
    #ax.set_zlim(-10., 10.)
    #ax.set_zticks(np.linspace(-10., 10., 5))
    #ax.set_zlabel(r'$z$', labelpad=25)
    ax.set_title('Custom')
    cb = fig.colorbar(p)
    cb.set_label('Weight (1e-6)', labelpad=50., rotation=270.)
    plt.tight_layout()

    # plotting the extended run
    ax = fig.add_subplot(122, projection='3d')
    p = ax.scatter(complete_results.samples[:, 0], complete_results.samples[:, 1], complete_results.samples[:, 2],
                   marker='o', c=np.exp(complete_results.logwt) * 1e8, linewidths=(0.,), cmap='coolwarm')
    #ax.set_xlim(-10., 10.)
    #ax.set_xticks(np.linspace(-10., 10., 5))
    #ax.set_xlabel(r'$x$', labelpad=25)
    #ax.set_ylim(-10., 10.)
    #ax.set_yticks(np.linspace(-10., 10., 5))
    #ax.set_ylabel(r'$y$', labelpad=25)
    #ax.set_zlim(-10., 10.)
    #ax.set_zticks(np.linspace(-10., 10., 5))
    #ax.set_zlabel(r'$z$', labelpad=25)
    ax.set_title('Complete')
    cb = fig.colorbar(p)
    cb.set_label('Weight (1e-8)', labelpad=50., rotation=270.)
    plt.tight_layout()

    plt.show()



def summary_plots(custom_results,
                  complete_results):

    # plot custom run
    fig, axes = dyplot.runplot(custom_results, color='red')

    # plot complete run
    fig, axes = dyplot.runplot(complete_results, color='dodgerblue', fig=(fig, axes))

    fig.tight_layout()

    plt.show()

with open('custom_sampler.dill','rb') as f:
    custom_sampler=dill.load(f)
    custom_results=custom_sampler.results
with open('complete_sampler.dill','rb') as f:
    complete_sampler=dill.load(f)
    complete_results=complete_sampler.results


summary_plots(custom_results,complete_results)
