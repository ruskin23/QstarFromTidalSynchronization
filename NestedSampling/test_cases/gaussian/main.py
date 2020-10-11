#Samples the 7 parameters of the main sampling
#We can impose a gaussian distribution on any of these


import argparse
import os
import sys

import matplotlib.pyplot as plt

import numpy
import scipy
from numpy import linalg
from scipy import stats


from dynesty import DynamicNestedSampler
from dynesty.dynamicsampler import stopping_function, weight_function

from pathos.pools import ProcessPool

import random
import copy

import pickle
import dill

def cmdline_args():

    parser = argparse.ArgumentParser()

    parser.add_argument('-s',
                        action='store',
                        dest='status',
                        help='start or continue sampling'
                        )


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




class NestedSampling:

    def loglikelihood(self,x):

        parameter_set=dict()
        for i,v in enumerate(self.sampling_parameters):
           parameter_set[v[0]]=x[i]

        L=0
        for key in self.model_parameters:
             L=L-0.5*(((parameter_set[key]-self.model_parameters[key]['value'])/self.model_parameters[key]['sigma'])**2)-numpy.log(self.model_parameters[key]['sigma']*numpy.sqrt(2*numpy.pi))

        return L

    def prior_transform(self,u):

        x=numpy.array(u)

        for i,s in enumerate(self.sampling_parameters):

            if s[-1]=='Normal':
                mean=s[1]
                sigma=s[2]
                x[i]=scipy.stats.norm.ppf(u[i],loc=mean,scale=sigma)
            elif s[-1]=='Truncated_normal':
                mean=s[1]
                sigma=s[2]
                low=(s[3]-mean)/sigma
                high=(s[4]-mean)/sigma
                x[i]=scipy.stats.truncnorm.ppf(u[i],low,high,loc=mean,scale=sigma)
            elif s[-1]=='Uniform':
                x[i]=(s[2]-s[1])*u[i] + s[1]

        return x



    def internal_sampling(self):

        dsampler=DynamicNestedSampler(self.loglikelihood, self.prior_transform,
                                       self.ndim,pool=self.pool,queue_size=self.queue_size)

        dsampler.run_nested()
        with open('internal_sampler.dill','wb') as f:
            dill.dump(dsampler,f)



    def nobreak_intial(self):
        #runs the full intiial samples loop and
        #generate dill file initial_samples.dill with all the initial samples

        dsampler=DynamicNestedSampler(self.loglikelihood, self.prior_transform,
                                       self.ndim,pool=self.pool,queue_size=self.queue_size)

        print_func=None
        print_progress=True
        pbar,print_func=dsampler._get_print_func(print_func,print_progress)
        niter=1
        ncall=0

        for results in dsampler.sample_initial():
            (worst, ustar, vstar, loglstar, logvol,
             logwt, logz, logzvar, h, nc, worst_it,
             boundidx, bounditer, eff, delta_logz) = results

            niter+=1
            ncall+=nc

            print_func(results,niter,ncall,nbatch=0,dlogz=0.01,logl_max=numpy.inf)

        with open('initial_samples.dill','wb') as f:
            dill.dump(dsampler,f)


    def break_initial(self,it,status):
        #breaks the initial sampler at iteration it
        #file saved as either initial_samples_*it*.dill or
        #intial_samples_end.dill

        if status=='start':
            dsampler=DynamicNestedSampler(self.loglikelihood, self.prior_transform,
                                       self.ndim,pool=self.pool,queue_size=self.queue_size)
            resume=False

        else:
            with open('initial_samples_100.dill','rb') as f:
                dsampler=dill.load(f)
            resume=True

            dsampler.rstate=numpy.random
            pool=ProcessPool(nodes=self.queue_size)
            dsampler.pool=pool
            dsampler.M=pool.map

            dsampler.sampler.rstate=numpy.random
            dsampler.sampler.pool=pool
            dsampler.sampler.M=pool.map

            dsampler.sampler.saved_h=dsampler.saved_h
            dsampler.sampler.saved_logz=dsampler.saved_logz
            dsampler.sampler.saved_logzvar=dsampler.saved_logzvar
            dsampler.sampler.saved_logvol=dsampler.saved_logvol
            dsampler.sampler.live_logl=dsampler.live_logl

        print_func=None
        print_progress=True
        pbar,print_func=dsampler._get_print_func(print_func,print_progress)

        niter=1
        ncall=0

        for results in dsampler.sample_initial(resume=resume):
            (worst, ustar, vstar, loglstar, logvol,
             logwt, logz, logzvar, h, nc, worst_it,
             boundidx, bounditer, eff, delta_logz) = results

            niter+=1
            ncall+=nc

            print_func(results,niter,ncall,nbatch=0,dlogz=0.01,logl_max=numpy.inf)

            if status=='start':
                if niter>it:
                    with open('initial_samples_'+str(it)+'.dill','wb') as f:
                        dill.dump(dsampler,f)
                    break

        if status=='continue':
            with open('initial_samples_end.dill','wb') as f:
                dill.dump(dsampler,f)


    def external_sampling(self,
                 status):


        self.nobreak_intial()
        #self.break_initial(100,status)

        #Adding batches to the initial sampling
        #dsampler.rstate=numpy.random
        #pool=ProcessPool(nodes=4)
        #dsampler.pool=pool
        #dsampler.M=pool.map

        #dsampler.sampler.rstate=numpy.random
        #dsampler.sampler.pool=pool
        #dsampler.sampler.M=pool.map


        #while True:
        #    stop = stopping_function(dsampler.results)  # evaluate stop
        #    if not stop:
        #        logl_bounds = weight_function(dsampler.results)  # derive bounds
        #        for results in dsampler.sample_batch(logl_bounds=logl_bounds):
        #            pass

        #        dsampler.combine_runs()  # add new samples to previous results
        #    else:
        #        break



    def __init__(self,
                 sampling_parameters,
                 pool,
                 queue_size):


        self.sampling_parameters=sampling_parameters
        self.pool=pool
        self.queue_size=queue_size
        self.ndim=7

        self.model_parameters=dict(logQ=dict(value=8.0,sigma=0.5),
                                   primary_mass=dict(value=0.8,sigma=0.1),
                                   age=dict(value=4.5,sigma=0.2))



args = cmdline_args()
system_number=args.system

status=args.status

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
            feh_value=float(x[4])
            feh_error=float(x[5])
            Porb_value=float(x[6])
            Porb_error=float(x[7])
            eccentricity_value=float(x[8])
            eccentricity_error=float(x[9])
            break



sampling_parameters = [('Porb',Porb_value,Porb_error,'Normal'),
                       ('feh',feh_value,feh_error,-1.014,0.537,'Truncated_normal'),
                       ('eccentricity',eccentricity_value,eccentricity_error,0.0,0.45,'Truncated_normal'),
                       ('Wdisk',2*scipy.pi/14,2*scipy.pi/1.4,'Uniform'),
                       ('logQ',5.0,12.0,'Uniform'),
                       ('primary_mass',0.5,1.2,'Uniform'),
                       ('age',1e-3,10.0,'Uniform')]



S=NestedSampling(sampling_parameters,pool,queue_size)
S.internal_sampling()
#S.external_sampling(status)
