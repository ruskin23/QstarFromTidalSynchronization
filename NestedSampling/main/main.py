import argparse
import os
import sys

poet_path='/home/ruskin/projects/poet/'
sys.path.append(poet_path+'PythonPackage')
sys.path.append(poet_path+'scripts')

import matplotlib.pyplot as plt

import numpy
from numpy import linalg
import scipy
from scipy import stats
import math
try:
    from scipy.special import logsumexp
except ImportError:
    from scipy.misc import logsumexp

import dynesty
from dynesty import plotting as dyplot
from dynesty.utils import *
from dynesty.dynamicsampler import *
from dynesty.results import *

from pathos.pools import ProcessPool

import random
import copy

import pickle
import dill

from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from stellar_evolution.derived_stellar_quantities import\
    TeffK,\
    LogGCGS,\
    RhoCGS
from evolution_class import evolution


current_directory=os.getcwd()
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




class DahiyaChutiya:

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



    def loglikelihood(self,x):

        #parameter_set=dict()
        #for i,v in enumerate(self.sampling_parameters):
        #   parameter_set[v[0]]=x[i]

        #model_set=self.calculate_model(x)
        #L=0
        #for key in self.observed_parameters:
        #    L=L-0.5*(((model_set[key]-self.observed_parameters[key]['value'])/self.observed_parameters[key]['sigma'])**2) - numpy.log(self.observed_parameters[key]['sigma']*numpy.sqrt(2*numpy.pi))


        #print('loglike: ',L)
        #return L

        # define likelihood constants
        ndim = 3
        C = numpy.identity(ndim)
        C[C==0] = 0.95
        Cinv = linalg.inv(C)
        lnorm = -0.5 * (numpy.log(2 * numpy.pi) * ndim + numpy.log(linalg.det(C)))

        y=numpy.array([x[4],x[5],x[6]])

        return -0.5 * numpy.dot(y, numpy.dot(Cinv, y)) + lnorm


    def prior_transform(self,u):

        x=numpy.array(u)

        for i,s in enumerate(self.sampling_parameters):

            if s[-1]=='Normal':
                mean=s[1]
                sigma=s[2]
                x[i]=scipy.stats.norm.ppf(u[i],loc=mean,scale=sigma)
            elif s[-1]=='Truncated_Normal':
                mean=s[1]
                sigma=s[2]
                low=(s[3]-mean)/sigma
                high=(s[4]-mean)/sigma
                x[i]=scipy.stats.truncnorm.ppf(u[i],low,high,loc=mean,scale=sigma)
            elif s[-1]=='Uniform':
                x[i]=(s[2]-s[1])*u[i] + s[1]

        return x


    def sample_initial(self,
                       dsampler,
                       print_func):

        live_points=None

        #Sample Initial Batch
        for results in dsampler.sample_initial(nlive=self.nlive_init,
                                               dlogz=self.dlogz_init,
                                               maxcall=self.maxcall_init,
                                               maxiter=self.maxiter_init,
                                               logl_max=self.logl_max_init,
                                               n_effective=self.n_effective_init,
                                               live_points=live_points):

            (worst, ustar, vstar, loglstar, logvol,
             logwt, logz, logzvar, h, nc, worst_it,
             boundidx, bounditer, eff, delta_logz) = results

            self.ncall+=nc
            self.niter+=1

            #print progress
            print_func(results,self.niter,self.ncall,nbatch=0,dlogz=self.dlogz_init,logl_max=self.logl_max_init)


            #with open(current_directory+'/dill_directory/system_84/sample_initial.dill','wb') as f:
            #    dill.dump(dsampler,f)
            #    dill.dump(self.ncall,f)
            #    dill.dump(self.niter,f)


        return results




    def add_batches(self,
                    dsampler,
                    print_func):


        for n in range(dsampler.batch,self.maxbatch):

            print('\nBatch = {}'.format(dsampler.batch))

            res=dsampler.results
            mcall=min(self.maxcall-self.ncall,self.maxcall_batch)
            miter=min(self.maxiter-self.niter,self.maxiter_batch)
            neff=dsampler.n_effective

            print('mcall = {}, miter = {}, neff = {}'.format(mcall,miter,neff))

            if mcall > 0  and miter > 0 and neff < numpy.inf:
                print('Evaluating Stopping Conditions')
                stop,stop_vals = stopping_function(res, self.stop_kwargs,
                                           rstate=dsampler.rstate, M=dsampler.M,
                                           return_vals=True)


                stop_post, stop_evid, stop_val = stop_vals
                print('Stopping Values: ',stop_vals)
            if mcall>0 and miter>0 and neff<numpy.inf and not stop:
                addBatch = dsampler.add_batch(nlive=self.nlive_batch,
                                              wt_function=self.wt_function,
                                              wt_kwargs=self.wt_kwargs,
                                              maxiter=miter,
                                              maxcall=mcall,
                                              save_bounds=self.save_bounds,
                                              print_progress=self.print_progress,
                                              print_func=print_func,
                                              stop_val=stop_val)


                self.ncall,self.niter,logl_bounds,results=addBatch
                with open(current_directory+'/dill_directory/system_84/sampler_at_batch_'+str(dsampler.batch)+'.dill','wb') as f:
                    dill.dump(dsampler,f)
                    dill.dump(self.ncall,f)
                    dill.dump(self.niter,f)
            elif logl_bounds[1]!=numpy.inf:
                break
            else:break



    def sampler(self):

        return dynesty.DynamicNestedSampler(self.loglikelihood, self.prior_transform,
                                       self.ndim,pool=self.pool,queue_size=self.queue_size)



    def __call__(self,
                 status):


        if status=='start':
            dsampler=self.sampler()
            print_func=None
            pbar,print_func=dsampler._get_print_func(print_func,self.print_progress)
            results=self.sample_initial(dsampler, print_func)

        elif status=='continue':
            with open(current_directory+'/dill_directory/system_84/sample_initial.dill','rb') as f:
                dsampler=dill.load(f)
                self.ncall=dill.load(f)
                self.niter=dill.load(f)

            dsampler.sampler.rstate=numpy.random
            dsampler.sampler.pool=self.pool
            dsampler.sampler.M=self.pool.map

            print_func=None
            pbar,print_func=dsampler._get_print_func(print_func,self.print_progress)

            self.add_batches(dsampler,
                             print_func)

        else: raise ValueError




    def __init__(self,
                 interpolator,
                 sampling_parameters,
                 observed_parameters,
                 fixed_parameters,
                 mass_ratio,
                 pool,
                 queue_size):


        self.interpolator=interpolator
        self.sampling_parameters=sampling_parameters
        self.observed_parameters=observed_parameters
        self.fixed_parameters=fixed_parameters
        self.mass_ratio=mass_ratio
        self.pool=pool
        self.queue_size=queue_size

        # Initialize values.
        self.maxcall = sys.maxsize
        self.maxiter = sys.maxsize
        self.maxiter_batch = sys.maxsize
        self.maxcall_batch = sys.maxsize
        self.maxbatch = sys.maxsize
        self.maxiter_init = sys.maxsize
        self.maxcall_init = sys.maxsize
        self.wt_function = weight_function
        self.wt_kwargs = dict()

        self.stop_function = None
        self.stop_kwargs = dict()
        self.save_bounds=True
        self.M=pool.map

        #initialize call values
        self.ndim=len(self.sampling_parameters)
        self.nlive_init=250
        self.nlive_batch=250
        self.dlogz_init=0.01
        self.n_effective_init=numpy.inf
        self.logl_max_init=numpy.inf
        self.ncall=0
        self.it=1
        self.niter=self.it-1
        self.logl_bounds=(-numpy.inf,numpy.inf)
        self.maxcall_init = min(self.maxcall_init, self.maxcall)  # set max calls
        self.maxiter_init = min(self.maxiter_init, self.maxiter)  # set max iterations

        self.print_progress=True


if __name__=='__main__':

    serialized_dir = poet_path +  "stellar_evolution_interpolators"
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    eccentricity_path=os.path.join(poet_path,'eccentricity_expansion_coef.txt').encode('ascii')

    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        eccentricity_path
    )


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



    sampling_parameters = [('Porb',Porb_value,Porb_error,'Normal'),
                           ('feh',feh_value,feh_error,-1.014,0.537,'Truncated_Normal'),
                           ('eccentricity',eccentricity_value,eccentricity_error,0.0,0.45,'Truncated_Normal'),
                           ('Wdisk',2*scipy.pi/14,2*scipy.pi/1.4,'Uniform'),
                           ('logQ',5.0,12.0,'Uniform'),
                           ('primary_mass',0.5,1.2,'Uniform'),
                           ('age',1e-3,10.0,'Uniform')]


    observed_parameters=dict(teff=dict(value=teff_value,
                                       sigma=teff_error),
                             logg=dict(value=logg_value,
                                       sigma=logg_error),
                             spin=dict(value=Pspin_value,
                                       sigma=Pspin_error)
                             )


    fixed_parameters = dict(disk_dissipation_age=5e-3,
                            planet_formation_age=5e-3,
                            wind=True,
                            wind_saturation_frequency=2.54,
                            diff_rot_coupling_timescale=5e-3,
                            wind_strength=0.17,
                            inclination=scipy.pi/2
                            )

    print('Sampling Parameters: ',sampling_parameters)
    print('Observed Parameters: ',observed_parameters)



    S=DahiyaChutiya(interpolator,
                    sampling_parameters,
                    observed_parameters,
                    fixed_parameters,
                    mass_ratio,
                    pool,
                    queue_size)
    S(status)
