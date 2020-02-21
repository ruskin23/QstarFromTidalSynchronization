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


    def _norm(self,
              value,
              loc=0.0,
              sigma=1.0):

        arg=(value-loc)/sigma
        return numpy.exp(-(arg**2)/2)

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


    def test_model(self,
                   parameter_set):


        L=1.0

        if self.test_case=='prior':return 1

        if self.test_case=='uncorrelated':
            for key in self.sampled_keys:
                param_value=parameter_set[key]
                L=L*self._norm(param_value,loc=self.sampled_parameters[key]['value'],sigma=self.model_width[key])
            return L

        if self.test_case=='correlated':
            for key in ['primary_mass','age']:
                L=L*self._norm(parameter_set[key],loc=self.sampled_parameters[key]['value'],sigma=self.model_width[key])
            arg=(parameter_set['primary_mass']*parameter_set['age']-self.sampling_parameters['primary_mass']['value']*self.sampling_parameters['age']['value'])/(self.model_width['primary_mass'])*(self.model_width['age'])
            return L*numpy.exp(-(arg**2)/2)

        if self.test_case=='correlated_logQ':
            for key in ['logQ','age']:
                L=L*self._norm(parameter_set[key],loc=self.sampling_parameters[key]['value'],sigma=self.model_width[key])
                print('L_key = ',L)
            L=L*self._norm(parameter_set['logQ']*parameter_set['age'],loc=self.sampling_parameters['logQ']['value']*self.sampling_parameters['age']['value'],sigma=self.model_width['logQ'])*(self.model_width['age'])

            print('L = ',L)
            return L


    def posterior_probability(self,
                              parameter_set=None):


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

        #Calculate Likelihood
        likelihood=self.test_model(parameter_set)
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
                S=S+float(x[3])*numpy.exp(-numpy.dot(arg,arg)/2)

        print('Transition Probability = ', S)
        posterior = (prior*likelihood)/S

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


    def propose_from_samples(self):

        N=0
        current_values=[self.current_parameters[key]/self.sampled_parameters[key]['step'] for key in
                        self.sampled_keys]

        print('\nCurrent Values = ', current_values)
        SampleSet=[]
        MM=[]

        with open(self.samples_file,'r') as f:
            next(f)
            for lines in f:
                x=lines.split()

                sample_set=[float(x[self.sampled_keys.index(key)]) for key in self.sampled_keys]
                sample_values=[float(x[self.sampled_keys.index(key)])/self.sampled_parameters[key]['step'] for key in
                        self.sampled_keys]


                mulitplicity=float(x[3])
                distance=numpy.subtract(sample_values,current_values)
                arg=numpy.dot(distance,distance)
                modified_mulitplicity=mulitplicity*numpy.exp(-arg)
                N=N+modified_mulitplicity

                MM.append(modified_mulitplicity)
                SampleSet.append(sample_set)

        U=random.uniform(0, 1)
        for s,mm in zip(SampleSet,MM):
            mm=mm/N
            if mm>U:return s
            else:U=U-mm

    def values_proposed(self):

        proposed = dict()

        for key,value in self.sampling_parameters.items():
            if value['dist']=='Normal' or value['dist']=='Uniform':
                proposed[key]=scipy.stats.norm.rvs(loc=self.current_parameters[key],scale=value['step'])


        s=self.propose_from_samples()

        proposed['primary_mass']=s[0]
        proposed['age']=s[1]
        proposed['feh']=s[2]

        for index,key in enumerate(self.sampled_keys):
            proposed[key]=s[index]


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
            else:f.write('\n')

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
                 output_directory,
                 test_case):

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
        self.sampled_parameters=dict()
        for key,value in self.sampling_parameters.items():
            if value['dist']=='Samples':
                self.sampled_parameters[key]=dict()
                while True:
                    self.sampled_parameters[key]['value']=value['value']
                    self.sampled_parameters[key]['step']=value['step']
                    break
                self.sampled_keys.append(key)
        print(self.sampled_parameters)
        self.model_width=dict(primary_mass=0.8,
                              age=0.8,
                              feh=0.2,
                              logQ=0.2)


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

        self.test_case=test_case

        self.time_stamp = 0.0
        self.time_elapsed = 0.0

