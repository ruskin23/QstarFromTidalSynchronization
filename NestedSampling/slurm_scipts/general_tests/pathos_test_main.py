import argparse
import os
import sys

import matplotlib.pyplot as plt

import numpy
from numpy import linalg
import scipy
from scipy import stats

import dynesty
from dynesty import plotting as dyplot

from pathos.pools import ProcessPool

import random


def cmdline_args():

    parser = argparse.ArgumentParser()

    parser.add_argument('-l',
                        action='store',
                        dest='system',
                        help='select a system for mcmc'
                        )

    parser.add_argument('-n',
                        action='store',
                        dest='num_threads',
                        help='number of parallel processes'
                        )


    return parser.parse_args()

class DahiyaChutiya:

    def loglikelihood(self,x):

        #sigma1=1.0
        #sigma2=0.5
        #rho=0.01

        #z=((x[6]-1.0)**2/sigma1) + ((x[4]-6.0)**2/sigma2) - (2*rho*(x[6]-1.0)*(x[4]-6.0)/sigma1*sigma2)

        #L = -(z/2*(1-rho**2)) - numpy.log(2*numpy.pi*sigma1*sigma2*numpy.sqrt(1-rho**2))

        # define likelihood constants
        ndim = 3
        C = numpy.identity(ndim)
        C[C==0] = 0.95
        Cinv = linalg.inv(C)
        lnorm = -0.5 * (numpy.log(2 * numpy.pi) * ndim + numpy.log(linalg.det(C)))

        y=numpy.array([x[4],x[5],x[6]])

        for i in range(1000000):
            u=random.random()

        return -0.5 * numpy.dot(y, numpy.dot(Cinv, y)) + lnorm


    def prior_transform(self,u):

        x=numpy.array(u)

        for i,s in enumerate(self.sampling_parameters):

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

    def sampler(self):

        ndim=len(self.sampling_parameters)

        dsampler=dynesty.NestedSampler(self.loglikelihood, self.prior_transform,
                                       ndim,nlive=1500,pool=self.pool,queue_size=queue_size)

        dsampler.run_nested()
        dresults=dsampler.results
        dresults.summary()

        # Plot a summary of the run.
        #rfig, raxes = dyplot.runplot(dresults)

        # Plot traces and 1-D marginalized posteriors.
        #tfig, taxes = dyplot.traceplot(dresults)

        # Plot the 2-D marginalized posteriors.
        #cfig, caxes = dyplot.cornerplot(dresults,show_titles=True)


        #rfig.tight_layout()
        #tfig.tight_layout()
        #cfig.tight_layout()
        #plt.show()



    def __init__(self,
                 sampling_parameters,
                 pool,
                 queue_size):

        self.sampling_parameters=sampling_parameters
        self.pool=pool
        self.queue_size=queue_size


args = cmdline_args()
system_number=args.system
num_threads=int(args.num_threads)
catalog_file='SpinlogQCatalog_el0.4.txt'

queue_size=num_threads
pool=ProcessPool(nodes=num_threads)

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


ndim=len(sampling_parameters)


S=DahiyaChutiya(sampling_parameters,pool,queue_size)
S.sampler()
