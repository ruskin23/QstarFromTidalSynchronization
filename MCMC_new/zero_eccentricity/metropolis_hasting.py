#"/home/kpenev/projects/git/poet/stellar_evolution_interpolators"

import scipy
from scipy.stats import norm
import numpy

import pickle
import csv

import os
import os.path


from evolution_class import evolution
from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library

import time
start_time = time.time()

class MetropolisHastings:

    def posterior_probability(self,parameter_set=None):
        """reutrns current surface spin of the star given the parameters"""

        find_spin = evolution(
                                self.interpolator,
                                parameter_set,
                                self.fixed_parameters,
                                self.mass_ratio,
                                self.instance
                               )

        self.spin_value = find_spin()

        if numpy.isnan(self.spin_value): return scipy.nan

        print('Current Spin Value = ', self.spin_value )
        prior = 1.0

        for (key_obs,value_obs),(key_parameter,value_parameter) in zip(self.observation_data.items(),parameter_set.items()):
            prior  *= scipy.stats.norm(value_obs['value'],value_obs['sigma']).pdf(value_parameter)
            print('prior = ', scipy.stats.norm(value_obs['value'],value_obs['sigma']).pdf(value_parameter))

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


        ext = 'age' + '\t' + 'current_orbital_period'  + '\t' + 'current_spin'  + '\t' + 'primary_mass'  + '\t' + 'secondary_mass' + '\n'

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
        load_ics_file = 'ics_data_' + self.instance + '.pickle'
        if os.path.isfile(load_solver_file) == True:
            with open(load_solver_file,'rb') as f:
                current_Porb=pickle.load(f)
                current_spin=pickle.load(f)
                primary_mass=pickle.load(f)
                secondary_mass=pickle.load(f)
                age=pickle.load(f)
        if os.path.isfile(load_ics_file)==True:
            with open(load_ics_file,'rb') as f1:
                binary_data=pickle.load(f1)
                initial_Porb=pickle.load(f1)
                current_Porb_discard=pickle.load(f1)
                current_spin_discard=pickle.load(f1)

            primary_envelop_angmom = binary_data['primary_envelope_angmom']
            primary_core_angmom = binary_data['primary_core_angmom']
            secondary_envelop_angmom = binary_data['secondary_envelope_angmom']
            secondary_core_angmom = binary_data['secondary_core_angmom']

        with open(f_name, 'a', 1) as file:
            file.write('%s\t' %self.iteration_step)
            for key, value in self.proposed_parameters.items():
                file.write('%s\t' %value)
            if os.path.isfile(load_solver_file) == True:
                file.write(
                    repr(age) + '\t' +
                    repr(current_Porb) + '\t' +
                    repr(current_spin) + '\t' +
                    repr(primary_mass) + '\t' +
                    repr(secondary_mass) + '\t'
                )
            if os.path.isfile(load_ics_file)==True:
                file.write(
                    repr(primary_envelop_angmom) + '\t' +
                    repr(primary_core_angmom) + '\t' +
                    repr(secondary_envelop_angmom) + '\t' +
                    repr(secondary_core_angmom) + '\t'
                )
            file.write( '\n')


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
        self.initial_parameters['Wdisk'] = numpy.random.uniform(low=self.Wdisk['min'],high=self.Wdisk['max'],size=None)
        self.initial_parameters['logQ'] = numpy.random.uniform(low=self.logQ['min'],high=self.logQ['max'],size=None)

        print ('\nINITIAL PARAMETERS SET')

        self.current_parameters = self.initial_parameters

        if self.iteration_step==1:self.write_header()

        print(self.current_parameters)



    def first_iteration(self):

        """Calculates first set of parameters and its corresponding posterior probability. The process will run until a non-zero or non-nan value of posterior probability is calculated"""

        while True:
            self.initialise_parameters()

            self.current_file_exist = None

            if self.current_parameters['feh']<-1.4014 or self.current_parameters['feh']>0.537:
                self.proposed_parameters = self.current_parameters
                self.isAccepted = False
                self.write_output()
                self.proposed_parameters = dict()
                self.iteration_step = self.iteration_step + 1
                print('feh value not in range')
                continue

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
            print('using accepted_paramter file')
            name = self.save_filename[0]
            self.current_file_exist = False
        else :
            print('using current_parameters file')
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
            print(step)

        self.iteration_step = int(array[0]) + 1

        self.current_parameters['teff_primary'] = float(array[1])
        self.current_parameters['feh'] = float(array[2])
        self.current_parameters['Porb'] = float(array[3])
        self.current_parameters['logg'] = float(array[4])
        self.current_parameters['Wdisk'] = float(array[5])
        self.current_parameters['logQ'] = float(array[6])

        if self.current_file_exist: self.current_posterior = float(array[7])
        else : self.current_posterior = self.posterior_probability(self.current_parameters)

        self.iterations()

    def __init__(
                self,
                interpolator,
                fixed_parameters,
                observation_data,
                Wdisk,
                logQ,
                proposed_step,
                total_iterations,
                observed_Pspin,
                mass_ratio,
                instance
                ):


        self.interpolator  = interpolator
        self.fixed_parameters = fixed_parameters
        self.observation_data = observation_data
        self.Wdisk=Wdisk
        self.logQ = logQ
        self.proposed_step = proposed_step
        self.iteration_step = 1
        self.total_iterations= total_iterations
        self.mass_ratio=mass_ratio

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
