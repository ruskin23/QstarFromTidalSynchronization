import argparse
import os
import sys

poet_path='/home/ruskin/projects/poet/'
sys.path.append(poet_path+'PythonPackage')
sys.path.append(poet_path+'scripts')

from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library

import numpy
import scipy
from scipy import stats
import matplotlib.pyplot as plt

import dynesty
from dynesty import utils as dyfunc
from dynesty import plotting as dyplot


def cmdline_args():

    parser = argparse.ArgumentParser()

    parser.add_argument('-l',
                        action='store',
                        dest='system',
                        help='select a system for mcmc'
                        )


    return parser.parse_args()

def loglike(x):

    #parameter_set=dict()
    #for i,key in enumerate(sampling_parameters):
    #    if key == 'age':print('age_i = ',i)
    #    if key == 'logQ':print('logQ_i = ',i)
    #    parameter_set[key]=x[i]

    L= (-0.5*(((x[6]-1.0)/0.5)**2)
        -0.5*(((x[4]-6.0)/0.2)**2)
        -numpy.log(0.5*numpy.sqrt(2*numpy.pi))
        -numpy.log(0.2*numpy.sqrt(2*numpy.pi))
        )
    return L


def ptform(u):

    x=numpy.array(u)

    for i,key in enumerate(sampling_parameters):

        if sampling_parameters[key]['dist']=='Normal':
            mean,sigma=sampling_parameters[key]['value'],sampling_parameters[key]['sigma']
            x[i]=scipy.stats.norm.ppf(u[i],loc=mean,scale=sigma)

        elif sampling_parameters[key]['dist']=='Truncated_Normal':
            mean,sigma=sampling_parameters[key]['value'],sampling_parameters[key]['sigma']
            low,high=(sampling_parameters[key]['low']-mean)/sigma,(sampling_parameters[key]['high']-mean)/sigma
            x[i]=scipy.stats.truncnorm.ppf(u[i],low,high,loc=mean,scale=sigma)

        elif sampling_parameters[key]['dist']=='Uniform':
            low,high=sampling_parameters[key]['low'],sampling_parameters[key]['high']
            x[i]=(high-low)*u[i]+low

        else:print('Insufficient Prior Information')


    return x







args = cmdline_args()
system_number=args.system

catalog_file='SpinlogQCatalog_el0.4.txt'

with open(catalog_file,'r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        at_system=x[0]
        if system_number==at_system:
            teff_value=float(x[2])
            teff_error=float(x[3])
            feh_value=float(x[4])
            feh_error=float(x[5])
            logg_value=float(x[10])
            logg_error=float(x[15])
            Porb_value=float(x[6])
            Porb_error=float(x[7])
            eccentricity_value=float(x[8])
            eccentricity_error=float(x[9])
            Pspin_value=float(x[12])
            Pspin_error=float(x[13])
            mass_ratio=float(x[14])
            break

sampling_parameters=dict(feh=dict(value=feh_value,
                                  sigma=feh_error,
                                  low=-1.014,
                                  high=0.537,
                                  dist='Truncated_Normal'),
                         Porb=dict(value=Porb_value,
                                   sigma=Porb_error,
                                   dist='Normal'),
                         eccentricity=dict(value=eccentricity_value,
                                           sigma=eccentricity_error,
                                           low=0.0,
                                           high=0.45,
                                           dist='Truncated_Normal'),
                         Wdisk=dict(low=2*scipy.pi/14,
                                    high=2*scipy.pi/1.4,
                                    dist='Uniform'),
                         logQ=dict(low=5.0,
                                   high=12.0,
                                   dist='Uniform'),
                         primary_mass=dict(low=0.5,
                                           high=1.2,
                                           dist='Uniform'),
                         age=dict(low=1e-3,
                                  high=10,
                                  dist='Uniform')

                         )



print('Sampling Parameters: ',sampling_parameters)

ndim=len(sampling_parameters)


print('\nStarting')
dsampler=dynesty.NestedSampler(loglike, ptform,
                               ndim,nlive=500)

rlist=[]

for i in range(3):
    dsampler.run_nested()
    rlist.append(dsampler.results)
    dsampler.reset()

# Merge into a single run.
dresults = dyfunc.merge_runs(rlist)

#print('Keys:', dresults.keys(),'\n')  # print accessible keys
#dresults.summary()  # print a summary
#print('\nRESULTS = ',dresults)

# Plot a summary of the run.
#rfig, raxes = dyplot.runplot(dresults)

# Plot traces and 1-D marginalized posteriors.
#tfig, taxes = dyplot.traceplot(dresults)

# Plot the 2-D marginalized posteriors.
cfig, caxes = dyplot.cornerplot(dresults,show_titles=True)


cfig.tight_layout()
plt.show()

