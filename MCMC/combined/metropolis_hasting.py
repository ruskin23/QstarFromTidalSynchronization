#!/usr/bin/env python3 -u

import scipy
from scipy.stats import norm
import numpy
import random

import pickle
import csv

import sys
import os
import os.path


from evolution_class import evolution
from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library

import time
start_time = time.time()

class MetropolisHastings:


    def write_header(self):

        f = self.output_directory
        file_list = os.listdir(f)
        for name in self.save_filename:
            if name in file_list:
                os.remove(name)

        header = []
        header.append('Iteration_step')

        for key in self.current_parameters.keys():
            header.append(key)

        ext = [ 'secondary_mass',
                'posterior_probability',
                'acceptance_ratio',
                'inital_orbital_period',
                'initial_eccentricity',
                'current_orbital_period',
                'current_eccentricity',
                'current_spin',
                'delta_e',
                'delta_p',
                'bad_solution_flag',
                'gsl_flag'
               ]


        for a in ext:
            header.append(a)

        for f_name in self.save_filename:
            with open(f_name, 'w', 1) as fheader:
                for name in header:
                    fheader.write('%s\t' % name)
                fheader.write('\n')

        with open(self.output_directory+'time_stamp_' + self.instance +'.txt','w',1) as f:
            f.write('Step'+'\t'+'Time_elapsed'+'\n')



    def posterior_probability(self,
                              parameter_set=None):


        model_calculations = evolution(self.interpolator,
                                       parameter_set,
                                       self.fixed_parameters,
                                       self.mass_ratio,
                                       self.instance,
                                       self.output_directory
                                       )

        self.spin=model_calculations()
        if numpy.isnan(self.spin):return scipy.nan

        print('Current Spin Value = ', self.spin)
        sys.stdout.flush()

        prior = 1.0

        #Calclate Priors
        for key,value in self.sampling_parameters.items():
            if value['dist']=='Normal':
                prior=prior*scipy.stats.norm(value['value'],value['sigma']).pdf(parameter_set[key])
                print('prior for {} = {}'.format(key,scipy.stats.norm(value['value'],value['sigma']).pdf(parameter_set[key])))
            if value['dist']=='Uniform':
                if numpy.logical_or(parameter_set[key]<self.sampling_parameters[key]['min'],
                                    parameter_set[key]>self.sampling_parameters[key]['max']):
                    prior=0
                    break

        #Calculate Likelihood
        print('observed Spin = ',self.observed_spin)
        likelihood=scipy.stats.norm(self.observed_spin['value'],self.observed_spin['sigma']).pdf(self.spin)
        print('likelihood = ', likelihood)
        sys.stdout.flush()

        #Calculate Transition Probability
        S=0
        phi_N=[parameter_set[key]/self.sampling_parameters[key]['step'] for key in self.sampled_keys]
        with open(self.samples_file,'r') as f:
            next(f)
            for lines in f:
                x=lines.split()
                phi_i=[float(x[self.sampled_keys.index(key)])/self.sampling_parameters[key]['step'] for key in self.sampled_keys]
                arg=numpy.subtract(phi_N,phi_i)
                S=S+numpy.exp(-float(x[3])*numpy.dot(arg,arg))

        print('Transition Probability = ', S)
        posterior = prior*likelihood*S

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
            sys.stdout.flush()


    def propose_from_sample(self,key):

        U=random.uniform(0, 1)
        parameter_value=self.current_parameters[key]
        with open(self.samples_file,'r') as f:
            next(f)
            for lines in f:
                x=lines.split()
                sample_value=float(x[self.sampled_keys.index(key)])
                mulitplicity=float(x[3])

                distance=parameter_value-sample_value
                modified_mulitplicity=mulitplicity*numpy.exp(-(distance/
                                                       self.sampling_parameters[key]['step'])**2)
                modified_parameter=sample_value*modified_mulitplicity
                if modified_parameter>U:return sample_value
                else:
                    U=U-modified_mulitplicity

    def values_proposed(self):

        proposed = dict()

        for key,value in self.sampling_parameters.items():
            if value['dist']=='Normal' or value['dist']=='Uniform':
                proposed[key]=scipy.stats.norm.rvs(loc=self.current_parameters[key],scale=value['step'])

            else:
                proposed[key]=self.propose_from_sample(key)

        self.proposed_parameters=proposed


    def write_output(self):

        if self.isAccepted == True:
            print('ACCEPTED')
            sys.stdout.flush()
            if self.iteration_step>1:
                self.current_parameters = self.proposed_parameters
                self.current_posterior = self.proposed_posterior
            f_name = self.save_filename[0]


        else:
            print('REJECTED')
            sys.stdout.flush()
            f_name = self.save_filename[1]


        load_solver_file = self.output_directory+'solver_results_'+self.instance+'.pickle'


        if os.path.isfile(load_solver_file) == True:
            with open(load_solver_file,'rb') as f:
                initial_orbital_period=pickle.load(f)
                intial_eccentricity=pickle.load(f)
                current_Porb=pickle.load(f)
                current_e=pickle.load(f)
                current_spin=pickle.load(f)
                delta_e=pickle.load(f)
                delta_p=pickle.load(f)
                bad_solution=pickle.load(f)
                gsl_flag=pickle.load(f)

        with open(f_name, 'a', 1) as f:

            f.write('%s\t' %self.iteration_step)

            if self.iteration_step==1:
                for key, value in self.current_parameters.items():
                    f.write('%s\t' % value)
                f.write(repr(self.current_parameters['primary_mass']*self.mass_ratio)+'\t'+
                        repr(self.current_posterior)+'\t'+
                        repr(self.p_acceptance)+'\t')

            else:
                for key, value in self.proposed_parameters.items():
                    f.write('%s\t' % value)
                f.write(repr(self.proposed_parameters['primary_mass']*self.mass_ratio)+'\t'+
                        repr(self.proposed_posterior)+'\t'+
                        repr(self.p_acceptance)+'\t')

            if os.path.isfile(load_solver_file) == True:
                f.write(
                    repr(initial_orbital_period)+'\t'+
                    repr(intial_eccentricity)+'\t'+
                    repr(current_Porb)+'\t'+
                    repr(current_e)+'\t'+
                    repr(current_spin)+'\t'+
                    repr(delta_e)+'\t'+
                    repr(delta_p)+'\t'+
                    repr(bad_solution)+'\t'+
                    repr(gsl_flag)+
                    '\n')

        with open(self.output_directory+'time_stamp_' + self.instance +'.txt', 'a',1) as f:
                f.write(repr(self.iteration_step) + '\t' + repr(self.time_elapsed) + '\n')

    def save_current_parameter(self):

        name = self.current_filename
        with open(name, 'w',1) as f:
            f.write(repr(self.iteration_step) + '\t')
            for key, value in self.current_parameters.items():
                f.write('%s\t' % value)
            f.write(repr(self.proposed_parameters['primary_mass']*self.mass_ratio)+'\t')
            f.write(repr(self.current_posterior)+ '\t' + repr(self.spin_value) +'\n')
        f.close()


    def initialise_parameters(self):

        """Initial values are drawn randomly from the given observable data set"""

        initial_parameters=dict()

        for key, value in self.sampling_parameters.items():
            initial_parameters[key]=value['value']

        with  open(self.solution_file) as f:
            next(f)
            for lines in f:
                x=lines.split()
                at_system=x[0]
                if at_system==self.system:
                    initial_parameters['logQ']=float(x[1])

        print ('\nINITIAL PARAMETERS SET:')
        for key,value in initial_parameters.items():
            print('{} = {}'.format(key,value))
        sys.stdout.flush()

        self.current_parameters = initial_parameters

        if self.iteration_step==1:self.write_header()


    def first_iteration(self):

        while True:

            self.initialise_parameters()
            self.current_file_exist = None
            self.current_posterior =  self.posterior_probability(parameter_set=self.current_parameters)
            if self.current_posterior==0 or numpy.isnan(self.current_posterior):
                continue
            else:
                break
        self.isAccepted = True
       	self.write_output()
       	self.iteration_step=self.iteration_step+1


    def iterations(self):

        """runs the specified number of iteration for the mcmc"""

        while True:

            self.isAccepted = None

            #initialising the set of parameters
            if self.iteration_step == 1 :
                self.first_iteration()
                continue

            #draw a random value from proposal function. The values will be proposed again if a negative age is encountered
            #while self.check_age_neg is True: self.values_proposed()
            self.values_proposed()
            print ('\nPROPOSED VALUES')
            print (self.proposed_parameters)
            sys.stdout.flush()
            if self.proposed_parameters['age']<0 or self.proposed_parameters['age']>10:
                self.isAccepted=False
                self.write_output()
                self.iteration_step=self.iteration_step+1
                print('age out of range')
                sys.stdout.flush()
                continue
            if self.proposed_parameters['primary_mass']<0.4 or self.proposed_parameters['primary_mass']>1.2:
                self.isAccepted=False
                self.write_output()
                self.iteration_step=self.iteration_step+1
                print('Mass out of range')
                sys.stdout.flush()
                continue
            if self.proposed_parameters['feh']<-1.014 or self.proposed_parameters['feh']>0.537:
                self.isAccepted = False
                self.write_output()
                self.iteration_step = self.iteration_step + 1
                print('feh value out of range')
                sys.stdout.flush()
                continue
            if self.proposed_parameters['eccentricity']<0:
                self.isAccepted = False
                self.write_output()
                self.iteration_step = self.iteration_step + 1
                print('e<0')
                sys.stdout.flush()
                continue


            #calculating posterior probabilty for proposed values. a nan value for posterior will mean the mass calculations were out of range of the interpolator. New values will be proposed.
            self.proposed_posterior = self.posterior_probability(parameter_set=self.proposed_parameters)
            if numpy.isnan(self.proposed_posterior):
                self.isAccepted = False
                self.write_output()
                self.iteration_step = self.iteration_step + 1
                print('Cannot calclate model parameters')
                sys.stdout.flush()
                continue


            # calculate acceptance probablity and check if proposed parameter is accepted or not
            self.check_acceptance()
            print ('isAccepted = ', self.isAccepted)
            sys.stdout.flush()

            self.time_elapsed = time.time() - self.time_stamp
            self.time_stamp = time.time()

            self.write_output()
            self.save_current_parameter()

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
            sys.stdout.flush()
        with open(self.save_filename[1], 'r') as f:
            reader = csv.reader(f, dialect='excel-tab')
            for row in reader:
                step = row
        self.iteration_step = int(step[0]) + 1

        self.current_parameters['primary_mass'] = float(array[5])
        self.current_parameters['age'] = float(array[6])
        self.current_parameters['feh'] = float(array[7])
        self.current_parameters['Porb'] = float(array[1])
        self.current_parameters['eccentricity'] = float(array[2])
        self.current_parameters['Wdisk'] = float(array[3])
        self.current_parameters['logQ'] = float(array[4])

        if self.current_file_exist: self.current_posterior = float(array[9])
        else : self.current_posterior = self.posterior_probability(self.current_parameters)

        self.iterations()

    def __init__(self,
                 system_number,
                 interpolator,
                 sampling_parameters,
                 fixed_parameters,
                 observed_spin,
                 catalog_file,
                 solution_file,
                 samples_file,
                 mass_ratio,
                 instance,
                 output_directory):

        self.system=system_number
        self.interpolator=interpolator
        self.sampling_parameters=sampling_parameters
        self.fixed_parameters=fixed_parameters
        self.iteration_step=1
        self.mass_ratio=mass_ratio
        self.catalog_file=catalog_file
        self.solution_file=solution_file
        self.samples_file=samples_file

        self.sampled_keys=[]
        for key,value in self.sampling_parameters.items():
            if value['dist']=='Samples':self.sampled_keys.append(key)

        self.current_parameters= dict()
        self.proposed_parameters = dict()

        self.isAccepted = None
        self.observed_spin=observed_spin
        self.check_age_neg = None

        self.proposed_posterior = 0.0
        self.current_posterior = 0.0
        self.p_acceptance = 0.0
        self.spin_value = 0.0

        self.instance = instance

        self.output_directory=output_directory

        self.save_filename = [
            self.output_directory+'accepted_parameters_'+self.instance+'.txt',
            self.output_directory+'rejected_parameters_'+ self.instance +'.txt'
            ]


        self.current_filename = self.output_directory+'current_parameters_'+ self.instance +'.txt'

        self.time_stamp = 0.0
        self.time_elapsed = 0.0

