#!/usr/bin/env python3
import sys
sys.path.append('home/ruskin/projects/poet/PythonPackage')
sys.path.append('home/ruskin/projects/poet/scripts')

from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as \
    orbital_evolution_library
from orbital_evolution.binary import Binary
from orbital_evolution.transformations import phase_lag
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet

from initial_secondary_angmom import IntialSecondaryAngmom
from initial_condition_solver import  InitialConditionSolver
from create_objects import BinaryObjects

import numpy
import scipy

wsun = 0.24795522138
#wsun=1.0

class Evolution:

    def calculate_intial_conditions(self):
        
        SecondaryAngmom=IntialSecondaryAngmom(self.interpolator,self.parameters)
        print('Seconary Initial Angular Momentum = ',SecondaryAngmom)()

        binary_system=BinaryObjects(self.interpolator,self.parameters)

        primary=binary_system.create_star(self.primary_mass,dissipation=True)
        secondary=binary_system.create_star(self.secondary_mass,dissipation=True)


        FindIC=InitialConditionSolver(self.interpolator,
                                      self.parameters,
                                      secondary_angmom=SecondaryAngmom)

        FindIC(primary,secondary)

    

    def evolve_binary(self,
                      parameter=None,
                      parameter_value=None):
        """A binary evolution function which creates primary, secondary and binary 
        system and evloves the binary according to specified parameters. 
        Args are optional. parameter passed in args will be assigned to self 
        
        Args:
            - parameter:
                name of a parameter if not specified while initializaing class
            - parameter_value:
                value of parameter passed
        
        Returns:
            -binary:
                a binary object after the evolution"""
        
        if parameter is not None:
            setattr(self,parameter,parameter_value)

        SecondaryAngmom=IntialSecondaryAngmom(self.interpolator,self.parameters)()
        print('Seconary Initial Angular Momentum = ',SecondaryAngmom)

        binary_system=BinaryObjects(self.interpolator,self.parameters)

        primary=binary_system.create_star(self.primary_mass,dissipation=True)
        secondary=binary_system.create_star(self.secondary_mass,dissipation=True)
        binary=binary_system.create_binary_system(primary,secondary,secondary_angmom=SecondaryAngmom)

        if self.print_cfile==True:
            if self.breaks==True:
                create_c_code='debug/cfile_'+self.system+'_withbreaks.cpp'
            else:
                create_c_code='debug/cfile_'+self.system+'.cpp'

            binary.evolve(
                self.age,
                self.evolution_max_time_step,
                self.evolution_precision,
                None,
                create_c_code=create_c_code,
                eccentricity_expansion_fname=b"eccentricity_expansion_coef.txt")
        else:
            binary.evolve(self.age,
                          self.evolution_max_time_step,
                          self.evolution_precision,
                          None)

        return binary

    def __init__(self,
                 interpolator,
                 parameters):

        self.interpolator=interpolator
        self.parameters=parameters

        for item,value in parameters.items():
            setattr(self,item,value)

        self.convective_phase_lag=phase_lag(self.logQ)
        print('Convective Phase lag = ',self.convective_phase_lag)



