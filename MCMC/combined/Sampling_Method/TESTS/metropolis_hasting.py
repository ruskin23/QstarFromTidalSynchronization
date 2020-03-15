#!/usr/bin/env python3 -u

import scipy
from scipy.stats import norm
import numpy
import random

import csv

import sys

sys.path.append('../Uncorrelated_Phi/')
from uncorrelated_sampling import UncorrelatedSampling

sys.path.append('../Adaptive/')
from adaptive_sampling import AdaptiveSampling
from covariance_matrix import Covariance
from utils import Norm
import os
import os.path

import time
start_time = time.time()

class MetropolisHastings:


    def test_model(self,
                   parameter_set):


        L=1.0

        normal=Norm()

        if self.test_case=='prior':return 1

        elif self.test_case=='gp':
            key=self.phi_key
            L=normal.N(parameter_set[self.phi_key],loc=self.sampled_parameters[key]['value'],sigma=self.model_width[key])
            return L

        else:
            if self.test_case in ['gpt','gptc']:key1=self.phi_key
            if self.test_case in ['gtt','gttc']:key1=self.theta_key
            key2='logQ'
            Norm1=normal.N(parameter_set[key1],loc=self.sampling_parameters[key1]['value'],sigma=self.model_width[key1])
            Norm2=normal.N(parameter_set[key2],loc=self.sampling_parameters[key2]['value'],sigma=self.model_width[key2])
            if self.test_case in ['gpt','gtt']:Norm3=1.0
            if self.test_case in ['gptc','gttc']:
                Norm3=numpy.exp(-(parameter_set[key1]-self.sampling_parameters[key1]['value'])*(parameter_set[key2]-self.sampling_parameters[key2]['value'])/(self.rho*self.model_width[key1]*self.model_width[key2]))
                #Norm3=numpy.exp(-(parameter_set[key1]*parameter_set[key2]-self.sampling_parameters[key1]['value']*self.sampling_parameters[key2]['value'])/(self.rho*self.model_width[key1]*self.model_width[key2]))
            L=Norm1*Norm2*Norm3
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

        #Calculate Stepping Probability
        phi_vector=numpy.array([parameter_set[key] for key in self.sampled_keys])
        if self.sampling_method=='uncorrelated':
            sampling=UncorrelatedSampling(self.system,
                                          self.sampling_parameters,
                                          parameter_set)
        if self.sampling_method=='adaptive':
            sampling=AdaptiveSampling(self.system,
                                      self.sampling_parameters,
                                      parameter_set,
                                      self.covariance_matrix)
        S=sampling.stepping_function_normalization(phi_vector)
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


    def propose_sample(self):

        proposed=dict()
        if self.sampling_method=='uncorrelated':
            proposed_samples=UncorrelatedSampling(self.system,self.sampling_parameters,self.current_parameters)
        if self.sampling_method=='adaptive':
            proposed_samples=AdaptiveSampling(self.system,self.sampling_parameters,self.current_parameters,self.covariance_matrix)
        proposed=proposed_samples()
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

        with open(f_name, 'a', 1) as f:

            f.write('%s\t' %self.iteration_step)

            if self.iteration_step==1:
                for key, value in self.current_parameters.items():
                    f.write('%s\t' % value)
                f.write(repr(self.current_posterior)+'\t'+
                        repr(self.p_acceptance)+'\n')

            else:
                for key, value in self.proposed_parameters.items():
                    f.write('%s\t' % value)
                f.write(repr(self.proposed_posterior)+'\t'+
                        repr(self.p_acceptance)+'\n')


    def save_current_parameter(self):

        name = self.current_filename
        with open(name, 'w',1) as f:
            f.write(repr(self.iteration_step) + '\t')
            for key, value in self.current_parameters.items():
                f.write('%s\t' % value)
            f.write(repr(self.current_posterior)+'\n')
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
        print(self.current_parameters)
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
            print('\nProposing Parameters from: ',self.current_parameters)
            self.propose_sample()
            print ('PROPOSED VALUES')
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

        for f_name in self.save_filename:
            with open(f_name, 'w', 1) as fheader:
                for name in header:
                    fheader.write('%s\t' % name)
                fheader.write('posterior'+'\t'+'p_acceptance'+'\n')




    def __init__(self,
                 system_number,
                 sampling_parameters,
                 instance,
                 sampling_method,
                 test_case,
                 output_directory,
                 catalog_file,
                 solution_file,
                 samples_file):

        self.system=system_number
        self.sampling_parameters=sampling_parameters
        self.instance=instance
        self.sampling_method=sampling_method
        self.test_case=test_case
        self.output_directory=output_directory
        self.catalog_file=catalog_file
        self.solution_file=solution_file
        self.samples_file=samples_file
        self.iteration_step=1

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

        self.phi_key='age'
        self.theta_key='eccentricity'
        self.model_width=dict(primary_mass=0.8,
                              age=0.8,
                              feh=0.2,
                              eccentricity=0.001,
                              logQ=0.2)

        self.rho=50

        self.current_parameters= dict()
        self.proposed_parameters = dict()

        self.isAccepted = None
        self.check_age_neg = None

        self.proposed_posterior = 0.0
        self.current_posterior = 0.0
        self.p_acceptance = 0.0

        self.save_filename = [
            self.output_directory+'accepted_parameters_'+self.instance+'.txt',
            self.output_directory+'rejected_parameters_'+ self.instance +'.txt'
            ]


        self.current_filename = self.output_directory+'current_parameters_'+ self.instance +'.txt'

        self.covariance_matrix=numpy.zeros([7,7])
        self.step_vector=numpy.zeros(7)
        k=0
        for key,item in self.sampling_parameters.items():
            self.step_vector[k]=item['step']
            k=k+1

        for i in range(7):
            for j in range(7):
                if i!=j:r=1/((i+1)*(j+1))
                else:r=1
                self.covariance_matrix[i,j]=r*self.step_vector[i]*self.step_vector[j]


        print('Covariance Matrix = ',self.covariance_matrix)



        #if sampling_method=='adaptive':
        #    self.covariance_matrix=Covariance(self.system).Calculate('Covariance')
