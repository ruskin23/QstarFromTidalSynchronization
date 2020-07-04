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

        #L= (-0.5*(((parameter_set['age']-1.0)/0.5)**2)
        #    -0.5*(((parameter_set['logQ']-6.0)/0.2)**2)
        #    -numpy.log(0.5*numpy.sqrt(2*numpy.pi))
        #    -numpy.log(0.2*numpy.sqrt(2*numpy.pi))
        #    )

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

        for i,s in enumerate(self.sampling_parameters):
            print('{} = {}'.format(s[0],x[i]))
        #for i,key in enumerate(self.sampling_parameters):

        #    if self.sampling_parameters[key]['dist']=='Normal':
        #        mean,sigma=self.sampling_parameters[key]['value'],self.sampling_parameters[key]['sigma']
        #        x[i]=scipy.stats.norm.ppf(u[i],loc=mean,scale=sigma)

        #    elif self.sampling_parameters[key]['dist']=='Truncated_Normal':
        #        mean,sigma=self.sampling_parameters[key]['value'],self.sampling_parameters[key]['sigma']
        #        low,high=(self.sampling_parameters[key]['low']-mean)/sigma,(self.sampling_parameters[key]['high']-mean)/sigma
        #        x[i]=scipy.stats.truncnorm.ppf(u[i],low,high,loc=mean,scale=sigma)

        #    elif self.sampling_parameters[key]['dist']=='Uniform':
        #        low,high=self.sampling_parameters[key]['low'],self.sampling_parameters[key]['high']
        #        x[i]=(high-low)*u[i]+low

        #    else:print('Insufficient Prior Information')


        return x


    def start(self):

        print('\nStarting')
        dsampler=dynesty.NestedSampler(self.loglike, self.ptform,
                                       self.ndim,nllive=500,pool=self.pool,queue_size=self.queue_size)
        dsampler.run_nested()
        dresults=dsampler.results
        dresults.summary()  # print a summary

        # Plot a summary of the run.
        rfig, raxes = dyplot.runplot(dresults)

        # Plot traces and 1-D marginalized posteriors.
        tfig, taxes = dyplot.traceplot(dresults)

        # Plot the 2-D marginalized posteriors.
        cfig, caxes = dyplot.cornerplot(dresults,show_titles=True)


        cfig.tight_layout()
        plt.savefig('cornerplot_84.png')



    def __init__(self,
                 system_number,
                 interpolator,
                 sampling_parameters,
                 fixed_parameters,
                 observed_parameters,
                 mass_ratio,
                 pool,
                 queue_size):
        self.system=system_number
        self.interpolator=interpolator
        self.sampling_parameters=sampling_parameters
        self.observed_parameters=observed_parameters
        self.fixed_parameters=fixed_parameters
        self.mass_ratio=mass_ratio

        self.ndim=len(self.sampling_parameters)
        self.pool=pool
        self.queue_size=queue_size


