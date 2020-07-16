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

class NestedSampling():


    def calculate_model(self,
                       x):

        parameter_set=dict()
        model_set=dict()

        for i,s in enumerate(self.sampling_parameters):
            parameter_set[s[0]]=x[i]

        print('\nParameter Set:')
        for key in parameter_set:
            print('{} = {}'.format(key,parameter_set[key]))

        spin_calculations = evolution(self.interpolator,
                                      parameter_set,
                                      self.fixed_parameters,
                                      self.mass_ratio)

        model_set['spin']=spin_calculations()

        mass=parameter_set['primary_mass']
        feh=parameter_set['feh']
        age=parameter_set['age']

        quantity_radius=self.interpolator('radius',mass,feh)
        quantity_lum=self.interpolator('lum',mass,feh)

        T=TeffK(quantity_radius,quantity_lum)
        try:model_set['teff']=T(parameter_set['age'])
        except:model_set['teff']=scipy.nan
        G=LogGCGS(mass,quantity_radius)
        try:model_set['logg']=G(parameter_set['age'])
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


        print('loglike: ',L)
        sys.stdout.flush()
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


    def start(self):

        print('\nStarting')
        sys.stdout.flush()
        dsampler=dynesty.DynamicNestedSampler(self.loglike, self.ptform, self.ndim,
                                              pool=self.pool, queue_size=self.queue_size)

        if self.sampling_type=='Internally':

            dsampler.run_nested()
            with open('complete_dynamic_sampler_results.dill','wb') as f:
                dill.dump(dsampler.results,f)

        else:

            #intialize sampler parameters
            nlive_init=100
            dlogz_init=0.01
            logl_max_init=numpy.inf
            n_effective_init=numpy.inf
            live_points=None

            print_func=None
            print_progress=True

            pbar,print_func=dsampler._get_print_func(print_func,print_progress)
            #Sample Initial Batch
            for results in dsampler.sample_initial(nlive=nlive_init,
                                                   dlogz=dlogz_init,
                                                   maxcall=self.maxcall_init,
                                                   maxiter=self.maxiter_init,
                                                   logl_max=logl_max_init,
                                                   n_effective=n_effective_init,
                                                   live_points=live_points):

                (worst, ustar, vstar, loglstar, logvol,
                 logwt, logz, logzvar, h, nc, worst_it,
                 boundidx, bounditer, eff, delta_logz) = results

                self.ncall+=nc
                self.niter+=1

                #print progress
                if print_progress:
                    print_func(results,self.niter,self.ncall,nbatch=0,dlogz=dlogz_init,logl_max=logl_max_init)

            res = dsampler.results
            with open(NestedSampling/slurm_scipts/system_tests+'/initial_sample_results.dill','wb') as f:
                dill.dump(dsampler.results,f)



    def __init__(self,
                 system_number,
                 interpolator,
                 sampling_parameters,
                 fixed_parameters,
                 observed_parameters,
                 mass_ratio,
                 pool,
                 queue_size,
                 sampling_type,
                 current_directory):

        self.system=system_number
        self.interpolator=interpolator
        self.sampling_parameters=sampling_parameters
        self.observed_parameters=observed_parameters
        self.fixed_parameters=fixed_parameters
        self.mass_ratio=mass_ratio
        self.sampling_type=sampling_type
        self.current_directory=current_directory
        # Initialize values.
        self.maxcall = sys.maxsize
        self.maxiter = sys.maxsize
        self.maxiter_batch = sys.maxsize
        self.maxcall_batch = sys.maxsize
        self.maxbatch = sys.maxsize
        self.maxiter_init = sys.maxsize
        self.maxcall_init = sys.maxsize
        self.wt_function = None
        self.wt_kwargs = dict()
        self.stop_function = None
        self.stop_kwargs = dict()

        #initialize call values
        self.ncall=0
        self.it=1
        self.niter=self.it-1
        self.logl_bounds=(-numpy.inf,numpy.inf)
        self.maxcall_init = min(self.maxcall_init, self.maxcall)  # set max calls

        self.ndim=len(self.sampling_parameters)
        self.pool=pool
        self.queue_size=queue_size



