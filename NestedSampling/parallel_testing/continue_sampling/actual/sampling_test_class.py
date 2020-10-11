import sys
import os
import numpy
import scipy
from scipy import stats

import dynesty
from dynesty import utils as dyfunc
from dynesty import plotting as dyplot
import matplotlib.pyplot as plt

from stellar_evolution.manager import StellarEvolutionManager
from stellar_evolution.derived_stellar_quantities import\
    TeffK,\
    LogGCGS,\
    RhoCGS

from evolution_class import evolution

import dill
import time

class NestedSampling():


    def calculate_model(self,
                       x):

        parameter_set=dict()
        model_set=dict()

        for i,s in enumerate(self.sampling_parameters):
            parameter_set[s[0]]=x[i]

        print('Parameter Set:')
        for key in parameter_set:
            print('{} = {}'.format(key,parameter_set[key]))

        mass=parameter_set['primary_mass']
        feh=parameter_set['feh']
        age=parameter_set['age']

        quantity_radius=self.interpolator('radius',mass,feh)
        quantity_lum=self.interpolator('lum',mass,feh)

        T=TeffK(quantity_radius,quantity_lum)
        try:model_set['teff']=T(age)
        except:model_set['teff']=scipy.nan

        G=LogGCGS(mass,quantity_radius)
        try:model_set['logg']=G(age)
        except:model_set['logg']=scipy.nan

        spin_calculations = evolution(self.interpolator,
                                      parameter_set,
                                      self.fixed_parameters,
                                      self.mass_ratio)

        model_set['spin']=spin_calculations()


        for key in model_set:
            if numpy.isnan(model_set[key]):model_set[key]=-scipy.inf
            print('{} = {}'.format(key,model_set[key]))

        return model_set


    def loglike(self,
                x):

        parameter_set=dict()
        for i,v in enumerate(self.sampling_parameters):
           parameter_set[v[0]]=x[i]

        model_set=self.calculate_model(x)
        L=0
        for key in self.observed_parameters:
            L=L-0.5*(((model_set[key]-self.observed_parameters[key]['value'])/self.observed_parameters[key]['sigma'])**2) - numpy.log(self.observed_parameters[key]['sigma']*numpy.sqrt(2*numpy.pi))


        print('loglike: ',L)
        return L

    def ptform(self,
               u):


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


    def SampleInitial(self,status):


        if status == 'start':
            dsampler=dynesty.DynamicNestedSampler(self.loglike, self.ptform, self.ndim,
                                              pool=self.pool, queue_size=self.queue_size)
            niter=1
            ncall=0
            resume=False

            #If initial sampling needs to be completely
            if self.sampling_type=='Internally':
                dsampler.run_nested()
                with open('complete_dynamic_sampler_results.dill','wb') as f:
                    dill.dump(dsampler.results,f)
                return


        elif status == 'continue':
            with open('initial_sampling_saved.dill','rb') as f:
                dsampler=dill.load(f)

            dsampler.loglikelihood=self.loglike
            dsampler.prior_transform=self.ptform

            dsampler.rstate=numpy.random
            dsampler.pool=self.pool
            dsampler.M=self.pool.map
            dsampler.queue_size=self.queue_size

            dsampler.sampler.loglikelihood=self.loglike
            dsampler.sampler.prior_transform=self.ptform

            dsampler.sampler.rstate=numpy.random
            dsampler.sampler.pool=self.pool
            dsampler.sampler.M=self.pool.map
            dsampler.sampler.queue_size=self.queue_size

            dsampler.sampler.saved_h=dsampler.saved_h
            dsampler.sampler.saved_logz=dsampler.saved_logz
            dsampler.sampler.saved_logzvar=dsampler.saved_logzvar
            dsampler.sampler.saved_logvol=dsampler.saved_logvol
            dsampler.sampler.live_logl=dsampler.live_logl


            niter=dsampler.it
            ncall=dsampler.ncall
            resume=True


        else:print('bhai sahi value dede')

        #Sample Initial Batch
        for results in dsampler.sample_initial(nlive=50,resume=resume):

            #print progress
            print_func=None
            print_progress=True
            pbar,print_func=dsampler._get_print_func(print_func,print_progress)

            (worst, ustar, vstar, loglstar, logvol,
             logwt, logz, logzvar, h, nc, worst_it,
             boundidx, bounditer, eff, delta_logz) = results

            ncall+=nc
            niter+=1

            print_func(results,niter,ncall,nbatch=0,dlogz=0.01,logl_max=numpy.inf)

            if status=='start':
                if niter>100:
                    print('Sampling paused {} iteration'.format(niter))
                    with open('initial_sampling_saved.dill','wb') as f:
                        dill.dump(dsampler,f)
                    return

        #save final results
        with open('initial_samples_finished.dill','wb') as f:
            dill.dump(dsampler,f)



    def __init__(self,
                 system_number,
                 interpolator,
                 sampling_parameters,
                 fixed_parameters,
                 observed_parameters,
                 mass_ratio,
                 pool,
                 queue_size,
                 sampling_type):

        self.system=system_number
        self.interpolator=interpolator
        self.sampling_parameters=sampling_parameters
        self.observed_parameters=observed_parameters
        self.fixed_parameters=fixed_parameters
        self.mass_ratio=mass_ratio
        self.sampling_type=sampling_type

        self.ndim=len(self.sampling_parameters)
        self.pool=pool
        self.queue_size=queue_size
        self.ndim=7


