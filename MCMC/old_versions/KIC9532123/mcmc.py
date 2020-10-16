#"/home/kpenev/projects/git/poet/stellar_evolution_interpolators"
import time
import pickle

import os
import os.path
import argparse
import csv

import sys
#sys.path.append('/Users/ruskinpatel/Desktop/Research/poet/PythonPackage')
#sys.path.append('/Users/ruskinpatel/Desktop/Research/poet/scripts')

sys.path.append('/home/kpenev/projects/git/poet/PythonPackage')
sys.path.append('/home/kpenev/projects/git/poet/scripts')

from evolution_class import evolution
from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library

import scipy
from scipy.stats import norm
import numpy

start_time = time.time()

class MetropolisHastings:

    def posterior_probability(self,parameter_set=None):
        """reutrns current surface spin of the star given the parameters"""

        find_spin = evolution(
                            self.interpolator,
                                parameter_set,
                                self.fixed_parameters,
                                self.instance
                               )

        self.spin_value = find_spin()

        if numpy.isnan(self.spin_value): return scipy.nan

        print('Current Spin Value = ', self.spin_value )
        prior = 1.0

        for (key_obs,value_obs),(key_parameter,value_parameter) in zip(self.observation_data.items(),parameter_set.items()):
            prior  *= scipy.stats.norm(value_obs['value'],value_obs['sigma']).pdf(value_parameter)

        likelihood = scipy.stats.norm(self.observed_Pspin['value'],self.observed_Pspin['sigma']).pdf(self.spin_value)


        print('likehood = ', likelihood)
        posterior = prior*likelihood
        return posterior


    def check_acceptance(self):
        self.p_acceptance = self.proposed_posterior/self.current_posterior
        if self.p_acceptance > 1: self.isAccepted = True
        else:
            rand = numpy.random.random_sample()
            if self.p_acceptance > rand : self.isAccepted = True
            else: self.isAccepted = False
            print('Acceptance_probability = ', self.p_acceptance)
            print('Random_number = ', rand)


    def values_proposed(self):

        proposed = dict()
        for (name_obs,value_obs),(name_step,value_step) in zip(self.current_parameters.items(),self.proposed_step.items()):
            proposed[name_obs]=scipy.stats.norm.rvs(loc=value_obs, scale=value_step)
            print("NAME AND VALUE",name_obs,proposed[name_obs] )

        self.proposed_parameters = proposed

    def write_header(self):

        f = os.getcwd()
        file_list = os.listdir(f)
        for name in self.save_filename:
            if name in file_list:
                os.remove(name)

        header = []
        header.append('Iteration_step')

        for key in self.current_parameters.keys():
            header.append(key)


        ext = 'inital_orbital_period'  + '\t' + 'initial_eccentricity'  + '\t' + 'current_orbital_period'  + '\t' + 'current_eccentricity'  + '\t' + 'current_spin'  + '\t' + 'delta_e'  + '\t' + 'delta_p'  + '\t' + 'primary_mass'  + '\t' + 'secondary_mass' + '\n'

        for f_name in self.save_filename:
            with open(f_name, 'w', 1) as file:
                for name in header:
                    file.write('%s\t' % name)
                file.write(ext)

    def write_output(self):


        if self.isAccepted == True:
            print('ACCEPTED')
            self.current_parameters = self.proposed_parameters
            self.current_posterior = self.proposed_posterior
            f_name = self.save_filename[0]


        else:
            print('REJECTED')
            f_name = self.save_filename[1]

        load_solver_file = 'solver_results_'+self.instance+'.pickle'

        if os.path.isfile(load_solver_file) == True:
            with open(load_solver_file,'rb') as f:
                ic=pickle.load(f)
                current_Porb=pickle.load(f)
                current_e=pickle.load(f)
                current_spin=pickle.load(f)
                delta_e=pickle.load(f)
                delta_p=pickle.load(f)
                primary_mass=pickle.load(f)
                secondary_mass=pickle.load(f)

        with open(f_name, 'a', 1) as file:
            file.write('%s\t' %self.iteration_step)
            for key, value in self.proposed_parameters.items():
                file.write('%s\t' % value)
            if os.path.isfile(load_solver_file) == True:
                file.write(
                    repr(ic[0]) + '\t' +
                    repr(ic[1]) + '\t' +
                    repr(current_Porb) + '\t' +
                    repr(current_e) + '\t' +
                    repr(current_spin) + '\t' +
                    repr(delta_e) + '\t' +
                    repr(delta_p) + '\t' +
                    repr(primary_mass) + '\t' +
                    repr(secondary_mass)  +
                    '\n')


    def save_current_parameter(self):

        name = self.current_filename
        with open(name, 'w') as f:
            f.write(repr(self.iteration_step) + '\t')
            for key, value in self.current_parameters.items():
                f.write('%s\t' % value)
            f.write(repr(self.current_posterior)+ '\t' + repr(self.spin_value) +'\n')
        f.close()


    def initialise_parameters(self):

        """Initial values are drawn randomly from the given observable data set"""


        for key, value in self.observation_data.items():
            self.initial_parameters[key] = scipy.stats.norm.rvs(loc=value['value'], scale=value['sigma'])
        self.initial_parameters['logQ'] = numpy.random.uniform(low=self.logQ['max'],high=self.logQ['max'],size=None)

        print ('\nINITIAL PARAMETERS SET')

        self.current_parameters = self.initial_parameters

        if self.iteration_step==1:self.write_header()

        print(self.current_parameters)



    def first_iteration(self):

        """Calculates first set of parameters and its corresponding posterior probability. The process will run until a non-zero or non-nan value of posterior probability is calculated"""

        while True:
            self.initialise_parameters()

            self.current_file_exist = None

            self.current_posterior =  self.posterior_probability(parameter_set=self.current_parameters)

            if self.current_posterior == 0 or numpy.isnan(self.current_posterior):
                self.proposed_parameters = self.current_parameters
                self.isAccepted = False
                self.write_output()
                self.proposed_parameters = dict()
                self.iteration_step = self.iteration_step + 1
                continue

            else: break


    def iterations(self):

        """runs the specified number of iteration for the mcmc"""


        #while self.iteration_step <= self.total_iterations:
        while True:

            self.isAccepted = None
            self.check_age_neg = False

            #initialising the set of parameters
            if self.iteration_step == 1 : self.first_iteration()

            #draw a random value from proposal function. The values will be proposed again if a negative age is encountered
            #while self.check_age_neg is True: self.values_proposed()
            self.values_proposed()
            print ('\nPROPOSED VALUES')
            print (self.proposed_parameters)
            if self.proposed_parameters['feh']<-1.014 or self.proposed_parameters['feh']>0.537:
                self.isAccepted = False
                self.write_output()
                self.iteration_step = self.iteration_step + 1
                print('feh value out of range')
                continue

            #calculating posterior probabilty for proposed values. a nan value for posterior will mean the mass calculations were out of range of the interpolator. New values will be proposed.
            self.proposed_posterior = self.posterior_probability(parameter_set=self.proposed_parameters)
            if numpy.isnan(self.proposed_posterior):
                self.isAccepted = False
                self.write_output()
                self.iteration_step = self.iteration_step + 1
                print("mass out of range for current values. proposing new values")
                continue

            # calculate acceptance probablity and check if proposed parameter is accepted or not
            print ('checking acceptance')
            self.check_acceptance()
            print ('isAccepted = ', self.isAccepted)


            self.time_elapsed = time.time() - self.time_stamp
            self.time_stamp = time.time()

            self.write_output()
            self.save_current_parameter()

            with open('time_stamp_' + self.instance +'.txt', 'a') as f:
                f.write(repr(self.iteration_step) + '\t' + repr(self.time_elapsed) + '\n')

            self.iteration_step = self.iteration_step + 1


    def continue_last(self):

        if  os.path.isfile(self.current_filename) == False or os.stat(self.current_filename).st_size == 0:
            name = self.save_filename[0]
            self.current_file_exist = False
        else :
            name = self.current_filename
            self.current_file_exist = True

        with open(name, 'r') as f:
            reader = csv.reader(f, dialect='excel-tab')
            for row in reader:
                array = row
            print(array)

        with open(self.save_filename[1], 'r') as f:
            reader = csv.reader(f, dialect='excel-tab')
            for row in reader:
                step = row
        self.iteration_step = int(step[0]) + 1

        self.current_parameters['teff_primary'] = float(array[1])
        self.current_parameters['feh'] = float(array[2])
        self.current_parameters['Porb'] = float(array[3])
        self.current_parameters['eccentricity'] = float(array[4])
        self.current_parameters['logg'] = float(array[5])
        self.current_parameters['Wdisk'] = float(array[6])
        self.current_parameters['logQ'] = float(array[7])

        if self.current_file_exist: self.current_posterior = float(array[8])
        else : self.current_posterior = self.posterior_probability(self.current_parameters)

        self.iterations()

    def __init__(
                self,
                interpolator,
                fixed_parameters,
                observation_data,
                logQ,
                proposed_step,
                total_iterations,
                observed_Pspin,
                instance
                ):


        self.interpolator  = interpolator
        self.fixed_parameters = fixed_parameters
        self.observation_data = observation_data
        self.logQ = logQ
        self.proposed_step = proposed_step
        self.iteration_step = 1
        self.total_iterations= total_iterations


        self.current_parameters= dict()
        self.proposed_parameters = dict()
        self.initial_parameters = dict()

        self.isAccepted = None
        self.observed_Pspin = observed_Pspin

        self.check_age_neg = None

        self.proposed_posterior = 0.0
        self.current_posterior = 0.0
        self.p_acceptance = 0.0
        self.spin_value = 0.0

        self.instance = instance

        #self.filename = ['accepted_test_1.txt', 'rejected_test_1.txt']
        self.save_filename = ['accepted_parameters_' + self.instance + '.txt',
                         'rejected_parameters_'+ self.instance +'.txt']


        self.current_filename = 'current_parameters_'+ self.instance +'.txt'

        self.time_stamp = 0.0
        self.time_elapsed = 0.0

if __name__ == '__main__':

    serialized_dir = "/home/ruskin/projects/poet/stellar_evolution_interpolators"
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        b"eccentricity_expansion_coef.txt"
    )


    observation_data = dict(
                        teff_primary=dict(value=5654, sigma=178),
                        feh=dict(value=-0.38, sigma=0.3),
                        Porb=dict(value=8.215, sigma=0.0001952),
                        eccentricity=dict(value=0.207, sigma=0.001),
                        logg=dict(value=4.6, sigma=0.14),
                        Wdisk=dict(value=2*scipy.pi / 1.4, sigma=0.1)
                    )

    observed_Pspin = dict(
                        value=7.197222223380878,
                        sigma=0.128 )

    fixed_parameters = dict(
                        disk_dissipation_age=5e-3,
                        planet_formation_age=5e-3,
                        wind=True,
                        wind_saturation_frequency=2.54,
                        diff_rot_coupling_timescale=5e-3,
                        wind_strength=0.17,
                        inclination=scipy.pi/2

    )

    proposed_step = dict(
                        teff_step=150.0,
                        feh_step=0.25,
                        Porb_step=0.0002,
                        eccentricity_step=0.001,
                        logg_step=0.1,
                        Wdisk_step=0.1,
                        logQ_step=0.08
                    )



    parser = argparse.ArgumentParser()
    parser.add_argument('-s', action='store_const', dest='start',
                    const='start',
                    help='start mcmc from beginning')

    parser.add_argument('-c', action='store_const', dest='cont',
                    const='continue',
                    help='continue mcmc from last iteration')

    parser.add_argument('-i', action = 'store', dest = 'instance',
                    help = 'define an instance of mcmc')
    parser.add_argument('-q',action = 'store', dest = 'logQ')
    args = parser.parse_args()

    logQ = dict(max=6.0 + int(args.instance)*0.083)


    #stepfilename = 'step_file_'+args.instance+'.txt'
    #with open(stepfilename,'w') as f:
        #for key,value in proposed_step.items():
        #f.write(key + '\t' + repr(value) + '\n')

    instance = args.instance
    mcmc = MetropolisHastings(
                            interpolator,
                            fixed_parameters,
                            observation_data,
                            logQ,
                            proposed_step,
                            10,
                            observed_Pspin,
                            instance)


    if args.start: mcmc.iterations()
    elif args.cont: mcmc.continue_last()
    else: print('provide correct arguments')
#flush()
#buffer_size =0 to write at the moment



#logQ = 0.05
#3sigma
#randomly select from masses

