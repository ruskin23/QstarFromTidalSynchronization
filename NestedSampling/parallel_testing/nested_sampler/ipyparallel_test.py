import argparse
import os
import sys

import scipy
import matplotlib.pyplot as plt

import ipyparallel as ipp

class Pool(object):

    def __init__(self, dview,nprocs):
        self.dview = dview
        self.size = nprocs

    def map(self, function, tasks):
        return self.dview.map_sync(function, tasks)

def cmdline_args():

    parser = argparse.ArgumentParser()

    parser.add_argument('-l',
                        action='store',
                        dest='system',
                        help='select a system for mcmc'
                        )


    return parser.parse_args()

def loglikelihood(x):

    sigma1=1.0
    sigma2=0.5
    rho=0.8

    z=((x[6]-1.0)**2/sigma1) + ((x[4]-6.0)**2/sigma2) - (2*rho*(x[6]-1.0)*(x[4]-6.0)/sigma1*sigma2)

    L = -(z/2*(1-rho**2)) - numpy.log(2*numpy.pi*sigma1*sigma2*numpy.sqrt(1-rho**2))

    return L


def prior_transform(u):

    x=numpy.array(u)

#    x[0]=scipy.stats.truncnorm.ppf(u[0],(sampling_parameters['feh']['low']-sampling_parameters['feh']['value'])/(sampling_parameters['feh']['sigma']),(sampling_parameters['feh']['high']-sampling_parameters['feh']['value'])/(sampling_parameters['feh']['sigma']),loc=sampling_parameters['feh']['value'],scale=sampling_parameters['feh']['sigma'])
#    x[1]=scipy.stats.norm.ppf(u[1],loc=sampling_parameters['Porb']['value'],scale=sampling_parameters['Porb']['sigma'])
#    x[2]=scipy.stats.truncnorm.ppf(u[2],(sampling_parameters['eccentricity']['low']-sampling_parameters['eccentricity']['value'])/(sampling_parameters['eccentricity']['sigma']),(sampling_parameters['eccentricity']['high']-sampling_parameters['eccentricity']['value'])/(sampling_parameters['eccentricity']['sigma']),loc=sampling_parameters['eccentricity']['value'],scale=sampling_parameters['eccentricity']['sigma'])
#    x[3]=(sampling_parameters['Wdisk']['high']-sampling_parameters['Wdisk']['low'])*u[3] + sampling_parameters['Wdisk']['low']
#    x[4]=(sampling_parameters['logQ']['high']-sampling_parameters['logQ']['low'])*u[4] + sampling_parameters['logQ']['low']
#    x[5]=(sampling_parameters['primary_mass']['high']-sampling_parameters['primary_mass']['low'])*u[5] + sampling_parameters['primary_mass']['low']
#    x[6]=(sampling_parameters['age']['high']-sampling_parameters['age']['low'])*u[6] + sampling_parameters['age']['low']



    for i,s in enumerate(sampling_parameters):

        if s[-1]=='Normal':
            mean=s[1]
            sigma=s[2]
            x[i]=scipy.stats.norm.ppf(u[i],loc=mean,scale=sigma)
        elif s[-1]=='Turncated_Normal':
            mean=s[1]
            sigma=s[2]
            low=(s[3]-mean)/sigma
            high=(s[4]-mean)/sigma
            x[i]=scipy.stats.truncnorm.ppf(u[i],low,high,loc=mean,scale=sigma)
        elif s[-1]=='Uniform':
            x[i]=(s[2]-s[1])*u[i] + s[1]



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
            break



sampling_parameters = [('Porb',Porb_value,Porb_error,'Normal'),
                       ('feh',feh_value,feh_error,-1.014,0.537,'Truncated_normal'),
                       ('eccentricity',eccentricity_value,eccentricity_error,0.0,0.45,'Truncated_normal'),
                       ('Wdisk',2*scipy.pi/14,2*scipy.pi/1.4,'Uniform'),
                       ('logQ',5.0,12.0,'Uniform'),
                       ('primary_mass',0.5,1.2,'Uniform'),
                       ('age',1e-3,10.0,'Uniform')]

#sampling_parameters=dict(feh=dict(value=feh_value,
#                                  sigma=feh_error,
#                                  low=-1.014,
#                                  high=0.537,
#                                  dist='Truncated_Normal'),
#                         Porb=dict(value=Porb_value,
#                                   sigma=Porb_error,
#                                   dist='Normal'),
#                         eccentricity=dict(value=eccentricity_value,
#                                           sigma=eccentricity_error,
#                                           low=0.0,
#                                           high=0.45,
#                                           dist='Truncated_Normal'),
#                         Wdisk=dict(low=2*scipy.pi/14,
#                                    high=2*scipy.pi/1.4,
#                                    dist='Uniform'),
#                         logQ=dict(low=5.0,
#                                   high=12.0,
#                                   dist='Uniform'),
#                         primary_mass=dict(low=0.5,
#                                           high=1.2,
#                                           dist='Uniform'),
#                         age=dict(low=1e-3,
#                                  high=10,
#                                  dist='Uniform')
#
#                         )

ndim=len(sampling_parameters)


#Starting Client
rc = ipp.Client()
nprocs = len(rc.ids)
print(rc.ids)

dview = rc[:]
dview.use_dill();

dview['sampling_parameters']=sampling_parameters

#sync imports on all engines
with dview.sync_imports():
    import numpy
    import scipy
    from scipy import stats
    #import random
    import dynesty
    from dynesty import plotting as dyplot

pool=Pool(dview,nprocs)

dsampler=dynesty.NestedSampler(loglikelihood, prior_transform,
                               ndim,nlive=5500,pool=pool)#,use_pool={'loglikelihood':False,'prior_transform':False})

dsampler.run_nested()
dresults=dsampler.results
dresults.summary()

# Plot the 2-D marginalized posteriors.
cfig, caxes = dyplot.cornerplot(dresults,show_titles=True)


cfig.tight_layout()
plt.show()

