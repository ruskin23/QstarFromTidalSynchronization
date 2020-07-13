
import argparse
import scipy
from scipy import stats
import numpy

import sys
import os
import os.path

from pathlib import Path
home_dir=str(Path.home())

if home_dir=='/home/rxp163130':
    poet_path=home_dir+'/poet/'
if home_dir=='/home/ruskin':
    poet_path=home_dir+'/projects/poet/'

sys.path.append(poet_path+'PythonPackage')


from stellar_evolution.manager import StellarEvolutionManager
from stellar_evolution.derived_stellar_quantities import\
    TeffK,\
    LogGCGS,\
    RhoCGS


class metropolis_hasting:

    def calculate_quantity(self,
                           name,
                           mass,
                           feh,
                           age):


        print('Calculating Quantity: ', name)
        quantity_radius=self.interpolator('radius',mass,feh)
        quantity_lum=self.interpolator('lum',mass,feh)

        if name=='teff':
            T=TeffK(quantity_radius,quantity_lum)
            print('Teff = ',T(age))
            return T(age)
        elif name=='logg':
            G=LogGCGS(mass,quantity_radius)
            print('Logg = ',G(age))
            return G(age)
        elif name=='rho':
            R=RhoCGS(mass,quantity_radius)
            return R(age)
        else:
            raise ValueError('quantity not found')


    def posterior_probability(self,
                              parameter_set):

        print('Before posterior CurrentResults = ', self.current_results)
        print('Calculating Posterior')

        evalute=dict()
        try:

            for key in self.constraints:
                evalute[key]=self.calculate_quantity(key,
                                                parameter_set['mass'],
                                                parameter_set['feh'],
                                                parameter_set['age'])
        except AssertionError:
            print('Maximum Age Error')
            return scipy.nan
        self.constraints_evaulated=evalute

        prior=scipy.stats.norm(self.parameters['feh']['value'],
                                 self.parameters['feh']['sigma']).pdf(parameter_set['feh'])
        print('Prior = ',prior)

        likelihood=1

        for key,val in self.constraints.items():
            likelihood=likelihood*scipy.stats.norm(self.constraints[key]['value'],
                                                   self.constraints[key]['sigma']).pdf(self.constraints_evaulated[key])
            print('likelihood = ', likelihood)

        posterior = likelihood*prior
        print('Posterior = ', posterior)
        print('after posterior CurrentResults = ', self.current_results)

        return posterior

    def check_acceptance(self):
        self.p_acceptance = self.proposed_posterior/self.current_posterior
        print('Acceptance Ratio = ', self.p_acceptance)
        print('in check_acceptance CurrentResults = ',self.current_results)
        if self.p_acceptance>1:
            print('Accepted')
            self.current_parameters=self.proposed_parameters
            self.current_posterior=self.proposed_posterior
            self.current_results=self.constraints_evaulated
            self.isAccepted=True
        else:
            print('Checking random number')
            rand=numpy.random.random_sample()
            print(rand)
            if self.p_acceptance>rand:
                print('Accepted')
                self.current_parameters=self.proposed_parameters
                self.current_posterior=self.proposed_posterior
                self.current_results=self.constraints_evaulated
            else:
                print('Rejected')
                self.isAccepted=False

    def propose_values(self):

        proposed=dict()
        print('proposing new parameters from:  ')
        print(self.current_parameters)
        for key in self.current_parameters:
            step_key=key+'_step'
            proposed[key]=scipy.stats.norm.rvs(loc=self.current_parameters[key],scale=self.step_size[step_key])

        self.proposed_parameters=proposed
        print('New parameters proposed:')
        print(self.proposed_parameters)

    def write_on_file(self,
                      fname,
                      parameter_set,
                      constraints_set):
        print('Writing On file')
        print('fname = ', fname)
        with open(fname,'a') as f:
            val=[repr(self.iteration_step)]
            for key,value in parameter_set.items():
                val.append(repr(value))
            for key,value in constraints_set.items():
                val.append(repr(value))
            f.write('\t'.join(val))
            f.write('\n')

    def write_output(self):

        if self.iteration_step==1:
            print('writing header')
            header=['Iteration_Step']
            for key in self.parameters:
                header.append(key)
            for key in self.constraints:
                header.append(key)
            for fname in self.filenames:
                with open(fname,'w+') as f:
                    f.write('\t'.join(header))
                    f.write('\n')

        print('in write output: CurrentResults = ', self.current_results)
        self.write_on_file(self.filenames[0],
                           self.current_parameters,
                           self.current_results)

        if self.isAccepted==False:
            self.write_on_file(self.filenames[1],
                               self.proposed_parameters,
                               self.constraints_evaulated)

    def initialise_parameters(self):

        for key,value in self.parameters.items():
            if key in ['mass','age']:
                self.current_parameters[key]=numpy.random.uniform(low=value['min'],high=value['max'])
            else:
                self.current_parameters[key]=scipy.stats.norm.rvs(loc=value['value'],scale=value['sigma'])


        if numpy.logical_or(self.current_parameters['feh']<-1.014,
                            self.current_parameters['feh']>0.537):
            print('Proposed Feh out of range')
            self.initialise_parameters()
        if numpy.logical_or(self.current_parameters['mass']<0.4,
                            self.current_parameters['mass']>1.2):
            print('Proposed mass out of range')
            self.initialise_parameters()
        if numpy.logical_or(self.current_parameters['age']<0.001,
                            self.current_parameters['age']>12.0):
            print('Proposed age out of range')
            self.initialise_parameters()

        print('First Parameter Set:')
        print(self.current_parameters)

        self.current_posterior=self.posterior_probability(self.current_parameters)

        self.current_results=self.constraints_evaulated

        if numpy.isnan(self.current_posterior):self.initialise_parameters()



    def iterations(self):
        #MCMC iterations

        #max_step=1000000

        while True:

            #Sample the first set of parameters and calculate its posterior
            if self.iteration_step==1:
                print('Initialising parameters')
                self.initialise_parameters()
                print('Current_posterior = ',self.current_posterior)
                self.write_output()
                self.iteration_step=self.iteration_step+1
                continue

            print('\n ITERATION STEP = ', self.iteration_step)

            #propose new parameters and calculate posterior
            self.propose_values()
            if numpy.logical_or(self.proposed_parameters['feh']<-1.014,
                                self.proposed_parameters['feh']>0.537):
                print('Proposed Feh out of range')
                self.isAccepted=False
                self.write_output()
                self.iteration_step=self.iteration_step+1
                continue
            if numpy.logical_or(self.proposed_parameters['mass']<0.4,
                                self.proposed_parameters['mass']>1.2):
                print('Proposed mass out of range')
                self.isAccepted=False
                self.write_output()
                self.iteration_step=self.iteration_step+1
                #if self.iteration_step>max_step:break
                continue
            if numpy.logical_or(self.proposed_parameters['age']<0.001,
                                self.proposed_parameters['age']>12.0):
                print('Proposed age out of range')
                self.isAccepted=False
                self.write_output()
                self.iteration_step=self.iteration_step+1
                #if self.iteration_step>max_step:break
                continue


            self.proposed_posterior=self.posterior_probability(self.proposed_parameters)
            print('Proposed_posterior = ', self.proposed_posterior)
            if numpy.isnan(self.proposed_posterior):
                isAccepted=False
                self.write_output()
                self.iteration_step=self.iteration_step+1
                #if self.iteration_step>max_step:break
                continue

            #calculate acceptance probability
            self.check_acceptance()

            #update output files
            self.write_output()

            #increment iteration step:
            self.iteration_step=self.iteration_step+1

            #if self.iteration_step>max_step:break



    def __init__(self,
                 interpolator,
                 parameters,
                 constraints,
                 step_size,
                 system):

        self.interpolator=interpolator
        self.parameters=parameters
        self.constraints=constraints
        self.step_size=step_size

        self.iteration_step=1

        self.current_parameters=dict()
        self.proposed_parameters=dict()
        self.constraints_evaulated=dict()
        self.current_results=dict()

        self.system=system

        self.isAccepted=None


        self.filenames=['samples/MassAgeFehSamples_'+self.system+'.txt',
                        'rejected_parameters/rejected_parameters_'+self.system+'.txt']

if __name__ == '__main__':

    parser=argparse.ArgumentParser()
    parser.add_argument('-l',action='store',dest='data_line',
                        help='specify the system from data_file.txt')
    args=parser.parse_args()

    system=args.data_line
    print('System = ', system)

    #interpolator
    serialized_dir = poet_path +  "stellar_evolution_interpolators"
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')


    with open('catalog_KIC.txt','r') as f:
        for lines in f:
            x=lines.split()
            at_system=x[0]
            if at_system==system:
                print('For values = ',lines)
                teff_value=float(x[2])
                teff_sigma=float(x[3])
                logg_value=float(x[10])
                logg_sigma=float(x[11])
                feh_value=float(x[4])
                feh_sigma=float(x[5])


    parameters = dict(
        mass=dict(min=0.5,max=1.2),
        age=dict(min=0.01,max=10),
        feh=dict(value=feh_value,sigma=feh_sigma)
    )

    print('Parameters = ',parameters)

    constraints = dict(
        teff=dict(value=teff_value,sigma=teff_sigma),
        logg=dict(value=logg_value,sigma=logg_sigma)
    )

    print('Constraints = ',constraints)

    step_size = dict(
        mass_step=1.0,
        age_step=1.0,
        feh_step=0.3
    )



    #initialise class
    MCMC=metropolis_hasting(interpolator,
                            parameters,
                            constraints,
                            step_size,
                            system)

    #metropolis_hasting
    MCMC.iterations()


