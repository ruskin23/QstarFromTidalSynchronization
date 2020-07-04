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


class NestedSampling:


    def calculate_model(self,
                       x):

        parameter_set=dict()
        model_set=dict()

        for i,key in enumerate(self.sampling_parameters):
            parameter_set[key]=x[i]

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

        #model_set=self.calculate_model(x)
        #L=0
        #for key in self.observed_parameters:
        #    L=L-0.5*(((model_set[key]-self.observed_parameters[key]['value'])/self.observed_parameters[key]['sigma'])**2) - numpy.log(self.observed_parameters[key]['sigma']*numpy.sqrt(2*numpy.pi))
        #print('lpglike = ',L)

        #TEST
        #age gaussian: model value = 1.0,width=0.5
        #logQ gaussian: model value = 6.0, width=0.2
        parameter_set=dict()
        for i,key in enumerate(self.sampling_parameters):
            parameter_set[key]=x[i]

        L= (-0.5*(((parameter_set['age']-1.0)/0.5)**2)
            -0.5*(((parameter_set['logQ']-6.0)/0.2)**2)
            -numpy.log(0.5*numpy.sqrt(2*numpy.pi))
            -numpy.log(0.2*numpy.sqrt(2*numpy.pi))
            )
        print('loglike = ',L)
        return L

    def ptform(self,
               u):


        x=numpy.array(u)

        for i,key in enumerate(self.sampling_parameters):

            if self.sampling_parameters[key]['dist']=='Normal':
                mean,sigma=self.sampling_parameters[key]['value'],self.sampling_parameters[key]['sigma']
                x[i]=scipy.stats.norm.ppf(u[i],loc=mean,scale=sigma)

            elif self.sampling_parameters[key]['dist']=='Truncated_Normal':
                mean,sigma=self.sampling_parameters[key]['value'],self.sampling_parameters[key]['sigma']
                low,high=(self.sampling_parameters[key]['low']-mean)/sigma,(self.sampling_parameters[key]['high']-mean)/sigma
                x[i]=scipy.stats.truncnorm.ppf(u[i],low,high,loc=mean,scale=sigma)

            elif self.sampling_parameters[key]['dist']=='Uniform':
                low,high=self.sampling_parameters[key]['low'],self.sampling_parameters[key]['high']
                x[i]=(high-low)*u[i]+low

            else:print('Insufficient Prior Information')

        #for i,key in enumerate(self.sampling_parameters):
        #    print('for parameter={}: u={} , ptform={}'.format(key,u[i],x[i]))


        return x

    def start(self):

        print('\nStarting')
        dsampler=dynesty.NestedSampler(self.loglike, self.ptform,
                                       self.ndim,nllive=1500)
        dsampler.run_nested(dlogz=0.1)
        dresults=dsampler.results
        print('Keys:', dresults.keys(),'\n')  # print accessible keys
        dresults.summary()  # print a summary
        #print('\nRESULTS = ',dresults)

        # Plot a summary of the run.
        rfig, raxes = dyplot.runplot(dresults)

        # Plot traces and 1-D marginalized posteriors.
        tfig, taxes = dyplot.traceplot(dresults)

        # Plot the 2-D marginalized posteriors.
        cfig, caxes = dyplot.cornerplot(dresults,show_titles=True)


        cfig.tight_layout()
        plt.show()



    def __init__(self,
                 system_number,
                 interpolator,
                 sampling_parameters,
                 fixed_parameters,
                 observed_parameters,
                 mass_ratio):

        self.system=system_number
        self.interpolator=interpolator
        self.sampling_parameters=sampling_parameters
        self.observed_parameters=observed_parameters
        self.fixed_parameters=fixed_parameters
        self.mass_ratio=mass_ratio

        self.ndim=len(self.sampling_parameters)
