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

from uncorrelated_sampling import UncorrelatedSampling


from evolution_class import evolution
from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library


class ModelTesting:

    def ModelCalculations(self,
                          parameter_set=None):

        parameter_array=[]
        for key in parameter_set:
            parameter_array=numpy.append(parameter_array,parameter_set[key])

        with open('model_16.pickle','rb') as f:
            s=pickle.load(f)
        parameter_array=numpy.append(parameter_array,parameter_set['primary_mass']*self.mass_ratio)
        result=s.predict(parameter_array.reshape(1,-1))
        self.Porb_iML,self.ecc_iML=result[0][0],result[0][1]

    def POETCalculations(self,
                    parameter_set=None):


        for key,value in self.sampling_parameters.items():
            if value['dist']=='Uniform':
                if numpy.logical_or(parameter_set[key]<self.sampling_parameters[key]['min'],
                                    parameter_set[key]>self.sampling_parameters[key]['max']):
                    return scipy.nan
                    break


        solve_for_spin = evolution(self.interpolator,
                                   parameter_set,
                                   self.fixed_parameters,
                                   self.mass_ratio,
                                   self.instance,
                                   self.output_directory
                                   )

        self.spin=solve_for_spin()
        if numpy.isnan(self.spin) or numpy.isinf(self.spin):return scipy.nan


    def propose_parameters(self):

        proposed = dict()

        self.proposing=UncorrelatedSampling(self.system,self.samples_file,self.sampling_parameters,self.proposed_parameters)
        proposed=self.proposing()
        self.proposed_parameters=proposed

    def write_header(self):

        f_name='system_output_'+self.instance+'.txt'
        header = []
        header.append('Iteration_step')

        for key in self.sampling_parameters.keys():
            header.append(key)

        ext = [ 'secondary_mass',
                'Porb_iPOET',
                'ecc_iPOET',
                'Porb_fPOET',
                'ecc_fPOET',
                'spinPOET',
                'delta_e',
                'delta_p',
                'bad_solution_flag',
                'gsl_flag',
                'binary_destroyed'
               ]


        for a in ext:
            header.append(a)

        with open(f_name, 'w', 1) as fheader:
            for name in header:
                fheader.write('%s\t' % name)
            fheader.write('\n')

    def write_output(self):


        f_name='system_output_'+self.instance+'.txt'

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

            for key, value in self.proposed_parameters.items():
                f.write('%s\t' % value)
            f.write(repr(self.proposed_parameters['primary_mass']*self.mass_ratio)+'\t')

            f.write(repr(self.Porb_iML)+'\t'+
                    repr(self.ecc_iML)+'\t')


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



    def iterations(self):

        while True:

            if self.iteration_step==1:
                proposed=dict()
                for key,value in self.sampling_parameters.items():
                    proposed[key]=value['value']
                self.proposed_parameters=proposed
                print('\nSTARTING PROPOSED VALUES')
                print (self.proposed_parameters)
                sys.stdout.flush()
                self.write_header()

            self.propose_parameters()
            print ('\nPROPOSED VALUES')
            print (self.proposed_parameters)
            sys.stdout.flush()

            self.POETCalculations(self.proposed_parameters)

            self.ModelCalculations(self.proposed_parameters)

            self.write_output()

            self.iteration_step = self.iteration_step + 1


    def __init__(self,
                 system_number,
                 interpolator,
                 samples_file,
                 sampling_parameters,
                 fixed_parameters,
                 mass_ratio,
                 instance,
                 current_directory,
                 output_directory):

        self.system=system_number
        self.interpolator=interpolator
        self.samples_file=samples_file
        self.sampling_parameters=sampling_parameters
        self.fixed_parameters=fixed_parameters
        self.iteration_step=1
        self.mass_ratio=mass_ratio
        self.Porb_iML=0.0
        self.ecc_iML=0.0


        self.phi_keys=[]
        self.theta_keys=[]
        for key,value in self.sampling_parameters.items():
            if value['dist'] in ['Normal','Uniform']:
                self.theta_keys.append(key)
            if value['dist']=='Samples':
                self.phi_keys.append(key)
        self.check_age_neg = None

        self.instance = instance

        self.output_directory=output_directory


