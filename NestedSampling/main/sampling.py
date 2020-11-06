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

from spin_calculation import SpinPeriod

import dill
import time
from datetime import date

import logging


today=date.today()
today=today.strftime("%B_%d")

class NestedSampling():


    def calculate_model(self,x,logger):

        parameter_set=dict()
        model_set=dict()

        for i,s in enumerate(self.sampling_parameters):
            parameter_set[s[0]]=x[i]

        logger.info('Parameter Set')
        for key in parameter_set:
            logger.info('{} = {}'.format(key,parameter_set[key]))

        spin_calculations = SpinPeriod(
                                    self.system,
                                    self.interpolator,
                                    parameter_set,
                                    self.fixed_parameters,
                                    self.mass_ratio,
                                    logger
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
            logger.info(f'{key} = {model_set[key]}')

        return model_set


    def loglike(self,
                x):

        logger=logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)
        pid=str(os.getpid())[-3:]
        scratch_filename=f'{self.output_directory}/system_{self.system}/{today}_sampling_{pid}.log'
        proces_handler=logging.FileHandler(scratch_filename)
        
        process_format=logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s',datefmt='%d-%b-%y %H:%M:%S')
        proces_handler.setFormatter(process_format)

        logger.addHandler(proces_handler)

        parameter_set=dict()
        for i,v in enumerate(self.sampling_parameters):
           parameter_set[v[0]]=x[i]

        model_set=self.calculate_model(x,logger)
        L=0
        for key in self.observed_parameters:
            L=L-0.5*(((model_set[key]-self.observed_parameters[key]['value'])/self.observed_parameters[key]['sigma'])**2) - numpy.log(self.observed_parameters[key]['sigma']*numpy.sqrt(2*numpy.pi))

        logger.info(f'Loglike = {L}')
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

    def initialize_sampler(self,dsampler):
        """If status=continue, assign attributes to the pickled sampler object"""

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

    def calculate_live_points(self,dsampler,nlive=500):
        """Samples nlive points from ndim unit cube, transform into desired priors and 
            calculates loglikelohoods of nlive points.
            uncomment section of code if the calulations is divided into various 
            instances and each intance is saved on a npz file if """
            
        #use pool to calculate points
        live_u=dsampler.rstate.rand(nlive,self.ndim)
        live_v=numpy.array(list(dsampler.M(self.ptform,numpy.array(live_u))))
        live_logl=numpy.array(list(dsampler.M(self.loglike,live_v)))
        
        #Convert -numpy.inf loglikelihoods to finite numbers
        for i,logl in enumerate(live_logl):
            if not numpy.isfinite(logl):
                if numpy.sign(logl)<0:
                    live_logl[i]=1e-300
                else:
                    raise ValueError(\
                        "The log-likelihood ({0}) of live "
                        "point {1} located at u={2} v={3} "
                        " is invalid."
                        .format(logl, i, live_u[i],live_v[i]))


        #If pool is not working or getting locked, use this method to create various
        #files and save results on npz files

        # live_logl=numpy.zeros(nlive)        
        # filename=self.output_directory+'/nlive_'+self.system+'.npz'
        # live_points=numpy.load(filename)

        # live_u=live_points['arr_0']
        # live_v=live_points['arr_1']
        # live_logl=live_points['arr_2']

        # for i,k in enumerate(live_logl):
        #     if k==0:break
        # rerun_idx=i
        # division=int((500-rerun_idx)/32)
        # instance=self.instance
        # i_min=rerun_idx + division*(instance-1)
        # if instance==32:i_max=500
        # else:i_max=rerun_idx + division*(instance)

        # print('Running between index {} and {}'.format(i_min,i_max))

        # for i in range(i_min,i_max):
        #     live_logl[i]=self.loglike(live_v[i])
        #     outfile=self.output_directory+'/nlive_'+self.system+'_'+str(instance)
        #     numpy.savez(outfile,live_u,live_v,live_logl)
            

    def get_live_points(self):
        """If live points were saved on a npz file, create a list of
        3 nd numpy arrays and conver -inf logl points tp finite numbers"""

        filename=self.results_directory+'/update_nlive/updated_nlive_'+self.system+'.npz'
        live_points=numpy.load(filename)

        live_u=live_points['arr_0']
        live_v=live_points['arr_1']
        live_logl=live_points['arr_2']

        #Convert -numpy.inf loglikelihoods to finite numbers
        for i,logl in enumerate(live_logl):
            if not numpy.isfinite(logl):
                if numpy.sign(logl)<0:
                    live_logl[i]=1e-300
                else:
                    raise ValueError(\
                        "The log-likelihood ({0}) of live "
                        "point {1} located at u={2} v={3} "
                        " is invalid."
                        .format(logl, i, live_u[i],live_v[i]))

        return [live_u,live_v,live_logl]

    def get_sampler_object(self,status):
        """Initialize Dynamic Nested Sampler object depending upon status"""
    
        if status == 'start':
            dsampler=dynesty.DynamicNestedSampler(self.loglike, self.ptform, self.ndim,pool=self.pool, queue_size=self.queue_size)
            self.niter=1
            self.ncall=0
            self.resume=False

        elif status == 'continue':
            with open(self.results_directory+'/initial_sampling_saved_'+str(self.system)+'.dill','rb') as f:
                dsampler=dill.load(f)

            dsampler=self.initialize_sampler(dsampler)
            self.niter=dsampler.it
            self.ncall=dsampler.ncall
            self.resume=True

        return dsampler

    def SampleInitial(self,status):


        dsampler=self.get_sampler_object(status)
        self.logger.info('Obtaining live points')
        live_points=self.get_live_points()
        self.logger.info('Live points obtained')
        #pylint: disable=unused-variable
        #Sample Initial Batch
        for results in dsampler.sample_initial(live_points=live_points,resume=self.resume):
            
            self.ncall+=results[9]
            self.niter+=1   

            result_names=['worst','ustar', 'vstar','loglstar','logvol','logwt', 'logz','logzvar','h','nc','worst_it','boundidx','bounditer', 'eff','delta_logz']
            self.logger.info(f'niter={self.niter}')
            self.logger.info(f'ncall={self.ncall}')
            for x,y in zip(result_names,results):
                self.logger.info(f'{x} = {y}')

            #overwrite. no need for time conditin. dump alternatively.
            with open(self.results_directory+'/initial_sampling_saved_'+str(self.system)+'.dill','wb') as f:
                dill.dump(dsampler,f)
            
            dsampler=self.initialize_sampler(dsampler)
        #pylint: enable=unused-variable

        #save final results
        with open(self.results_directory+'/initial_samples_finished_'+str(self.system)+'.dill','wb') as f:
            dill.dump(dsampler,f)



    def __init__(self,
                 system_number,
                 instance,
                 interpolator,
                 sampling_parameters,
                 fixed_parameters,
                 observed_parameters,
                 mass_ratio,
                 pool,
                 queue_size,
                 results_directory,
                 output_directory,
                 logger):

        self.system=system_number
        self.instance=instance
        self.interpolator=interpolator
        self.sampling_parameters=sampling_parameters
        self.observed_parameters=observed_parameters
        self.fixed_parameters=fixed_parameters
        self.mass_ratio=mass_ratio
        self.ndim=len(self.sampling_parameters)
        self.pool=pool
        self.queue_size=queue_size
        self.results_directory=results_directory
        self.output_directory=output_directory

        self.logger=logger
        


