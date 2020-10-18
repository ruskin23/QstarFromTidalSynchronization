#!/usr/bin/env python3 -u

import scipy
from scipy.stats import norm
import numpy
import random
import shelve

import pickle
import csv

import os
import os.path

import sys
from pathlib import Path
home_dir=str(Path.home())
git_dir='/QstarFromTidalSynchronization/MCMC/combined'
if home_dir=='/home/rxp163130':current_directory=home_dir+git_dir
if home_dir=='/home/ruskin':current_directory=home_dir+'/projects'+git_dir
if home_dir=='/home1/06850/rpatel23':
    work_dir='/work/06850/rpatel23/stampede2'
    current_directory=work_dir+git_dir
sys.path.append(current_directory+'/Sampling_Method/Uncorrelated_Phi')
from uncorrelated_sampling import UncorrelatedSampling
sys.path.append(current_directory+'/Sampling_Method/Adaptive')
from adaptive_sampling import AdaptiveSampling


from evolution_class import evolution
from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library

import time

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

        ext = [ 'posterior_probability',
                'acceptance_ratio',
                'inital_orbital_period',
                'initial_eccentricity',
                'current_orbital_period',
                'current_eccentricity',
                'current_spin',
                'delta_e',
                'delta_p',
                'bad_solution_flag',
                'gsl_flag',
                'binary_destroyed'
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


        prior = 1.0

        #Calclate Priors
        for key,value in self.sampling_parameters.items():
            if value['dist']=='Normal':
                prior=prior*scipy.stats.norm(value['value'],value['sigma']).pdf(parameter_set[key])
                print('prior for {} = {}'.format(key,scipy.stats.norm(value['value'],value['sigma']).pdf(parameter_set[key])))
            if parameter_set[key]<value['min_value'] or parameter_set[key]>value['max_value']:
                print('For parameter value {} = {}, proposed value out of range'.format(key,parameter_set[key]))
                prior=0
                break

        if prior==0: return scipy.nan

        model_calculations = evolution(self.interpolator,
                                       parameter_set,
                                       self.fixed_parameters,
                                       self.instance,
                                       self.output_directory
                                       )

        self.spin=model_calculations()
        if numpy.isnan(self.spin) or numpy.isinf(self.spin):return scipy.nan

        print('Current Spin Value = ', self.spin)
        sys.stdout.flush()


        #Calculate Likelihood
        print('observed Spin = ',self.observed_spin)
        likelihood=scipy.stats.norm(self.observed_spin['value'],self.observed_spin['sigma']).pdf(self.spin)
        print('likelihood = ', likelihood)
        sys.stdout.flush()

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
            sys.stdout.flush()


    def values_proposed(self):

        proposed = dict()

        for key,item in self.current_parameters.items():
            proposed[key]=scipy.stats.norm.rvs(loc=self.current_parameters[key],scale=self.sampling_parameters[key]['step'])

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
                binary_destroyed=pickle.load(f)

        with open(f_name, 'a', 1) as f:

            f.write('%s\t' %self.iteration_step)

            if self.iteration_step==1:posterior=self.current_posterior
            else:posterior=self.proposed_posterior
            for key, value in self.proposed_parameters.items():
                f.write('%s\t' % value)
            f.write(repr(posterior)+'\t'+
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
                    repr(gsl_flag)+'\t'+
                    repr(binary_destroyed)+
                    '\n')
            else:f.write('\n')

        with open(self.output_directory+'time_stamp_' + self.instance +'.txt', 'a',1) as f:
                f.write(repr(self.iteration_step) + '\t' + repr(self.time_elapsed) + '\n')

    def save_current_parameter(self):

        name = self.current_filename
        with open(name, 'w',1) as f:
            f.write(repr(self.iteration_step) + '\t')
            for key, value in self.sampling_parameters.items():
                f.write('%s\t' % self.current_parameters[key])
            f.write(repr(self.current_posterior)+ '\t' + repr(self.spin_value) +'\n')
        f.close()


    def initialise_parameters(self):

        """Initial values are drawn randomly from the given observable data set"""

        initial_parameters=dict()

        for key, value in self.sampling_parameters.items():
            initial_parameters[key]=value['value']

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

            self.time_stamp = time.time()

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

            self.write_output()
            self.save_current_parameter()

            self.iteration_step = self.iteration_step + 1


    def continue_last(self):

        with open(self.save_filename[1],'r') as f:
            line=list(f)[-1]
            x=line.split()
            iter_step_r=int(x[0])

        with open(self.save_filename[0],'r') as f:
            line=list(f)[-1]
            print(line)
            x=line.split()
            iter_step_a=int(x[0])
            self.current_parameters['Porb'] = float(x[1])
            self.current_parameters['ecosw'] = float(x[2])
            self.current_parameters['esinw'] = float(x[3])
            self.current_parameters['primary_mass'] = float(x[4])
            self.current_parameters['secondary_mass'] = float(x[5])
            self.current_parameters['age'] = float(x[6])
            self.current_parameters['feh'] = float(x[7])
            self.current_parameters['Wdisk'] = float(x[8])
            self.current_parameters['logQ'] = float(x[9])
            self.current_posterior = float(x[10])

        self.iteration_step=max(iter_step_r,iter_step_a) + 1

        print('\nContinuing from {} and {}'.format(self.iteration_step,self.current_parameters))

        self.iterations()


    def __init__(self,
                 system_number,
                 interpolator,
                 sampling_parameters,
                 fixed_parameters,
                 observed_spin,
                 instance,
                 current_directory,
                 output_directory):

        self.system=system_number
        self.interpolator=interpolator
        self.sampling_parameters=sampling_parameters
        self.fixed_parameters=fixed_parameters
        self.iteration_step=1

        self.current_parameters= dict()
        self.proposed_parameters = dict()

        self.isAccepted = None
        self.observed_spin=observed_spin

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

