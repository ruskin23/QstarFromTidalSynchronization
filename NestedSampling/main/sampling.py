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

#from evolution_class import evolution
from spin_calculation import SpinPeriod

import dill
import time

class NestedSampling:


    def calculate_model(self,
                       x):

        parameter_set=dict()
        model_set=dict()

        for i,s in enumerate(self.sampling_parameters):
            parameter_set[s[0]]=x[i]

        print('\nParameter Set:')
        for key in parameter_set:
            print('{} = {}'.format(key,parameter_set[key]))

        spin_calculations = SpinPeriod(
                                    self.system,
                                    self.interpolator,
                                    parameter_set,
                                    self.fixed_parameters,
                                    self.mass_ratio
                                    )

        model_set['spin']=spin_calculations()

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


        print('loglike: ',L, file=sys.stdout, flush=True)
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

        for i,s in enumerate(self.sampling_parameters):
            print('For {} = {}'.format(s[0],x[i]))

        return x

    def initialize_sampler(self,dsampler):

            dsampler.rstate=numpy.random
            dsampler.pool=self.pool
            dsampler.M=self.pool.map

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

            return dsampler


    def SampleInitial(self,status):

        with open('progress.txt','w',1) as f:
            f.write('niter'+'\t'+
                    'ncall'+'\t'+
                    'delta_logz'+'\t'+
                    'ustar'+'\t'+
                    'vstar'+'\t'+
                    'loglstar'+'\t'+
                    'logvol'+'\t'+
                    'logwt'+'\t'+
                    'logz'+'\t'+
                    'logzvar'+'\t'+
                    'worst_it'+'\t'+
                    'boundidx'+'\t'+
                    'bounditer'+'\t'+
                    'eff'+'\n')

        if status == 'start':
            dsampler=dynesty.DynamicNestedSampler(self.loglike, self.ptform, self.ndim,pool=self.pool, queue_size=self.queue_size)
            niter=1
            ncall=0
            resume=False

        elif status == 'continue':
            with open(self.output_directory+'/initial_sampling_saved.dill','rb') as f:
                dsampler=dill.load(f)

            dsampler=self.initialize_sampler(dsampler)
            niter=dsampler.it
            ncall=dsampler.ncall
            resume=True


        else:print('status can either be start or continue', file=sys.stdout, flush=True)

        #pylint: disable=unused-variable
        #Sample Initial Batch
        start_time=time.time()
        for results in dsampler.sample_initial(resume=resume):

            (worst, ustar, vstar, loglstar, logvol,
             logwt, logz, logzvar, h, nc, worst_it,
             boundidx, bounditer, eff, delta_logz) = results

            ncall+=nc
            niter+=1

            with open(self.output_directory+'/progress.txt','a',1) as f:
                f.write(repr(niter)+'\t'+
                        repr(ncall)+'\t'+
                        repr(delta_logz)+'\t'+
                        repr(ustar)+'\t'+
                        repr(vstar)+'\t'+
                        repr(loglstar)+'\t'+
                        repr(logvol)+'\t'+
                        repr(logwt)+'\t'+
                        repr(logz)+'\t'+
                        repr(logzvar)+'\t'+
                        repr(worst_it)+'\t'+
                        repr(boundidx)+'\t'+
                        repr(bounditer)+'\t'+
                        repr(eff)+'\n')

            #Stop sampling at 4 hours on cluster time
            time_spent=(time.time()-start_time)*0.000277778 #in hours
            print('Time spent = ',time_spent, file=sys.stdout, flush=True)

            if time_spent>4:
                print('Sampling paused {} iteration after spending {} hours'.format(niter,time_spent), file=sys.stdout, flush=True)
                with open(self.output_directory+'/initial_sampling_saved.dill','wb') as f:
                    dill.dump(dsampler,f)
                return
        #pylint: enable=unused-variable

        #save final results
        with open(self.output_directory+'/initial_samples_finished.dill','wb') as f:
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
                 output_directory):

        self.system=system_number
        self.interpolator=interpolator
        self.sampling_parameters=sampling_parameters
        self.observed_parameters=observed_parameters
        self.fixed_parameters=fixed_parameters
        self.mass_ratio=mass_ratio
        self.ndim=len(self.sampling_parameters)
        self.pool=pool
        self.queue_size=queue_size
        self.output_directory=output_directory


