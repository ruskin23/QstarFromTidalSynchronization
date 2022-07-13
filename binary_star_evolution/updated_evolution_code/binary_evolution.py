#!/usr/bin/env python3
from distutils.log import info
import sys
from pathlib import Path
from directories import directories
import logging


home_dir=str(Path.home())
path=directories(home_dir)
sys.path.append(path.poet_path+'/PythonPackage')
sys.path.append(path.poet_path+'/scripts')

from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as \
    orbital_evolution_library
from orbital_evolution.binary import Binary
from orbital_evolution.transformations import phase_lag
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet

from initial_secondary_angmom import IntialSecondaryAngmom
from initial_condition_solver import  InitialConditionSolver
#from intial_orbital_period_solver import InitialConditionSolver
#from initial_conditions_minimizer import InitialConditionSolver
from create_objects import BinaryObjects

import numpy
import scipy

wsun = 0.24795522138
#wsun=1.0

_logger=logging.getLogger(__name__)

class Evolution:

    def calculate_intial_conditions(self):

        SecondaryAngmom=IntialSecondaryAngmom(self.interpolator,self.parameters)
        SA=SecondaryAngmom()
        _logger.info('Seconary Initial Angular Momentum = {!r}'.format(repr(SA)))

        binary_system=BinaryObjects(self.interpolator,self.parameters)

        primary=binary_system.create_star(self.primary_mass,dissipation=True)
        secondary=binary_system.create_star(self.secondary_mass,dissipation=True)


        FindIC=InitialConditionSolver(self.interpolator,
                                      self.parameters,
                                      secondary_angmom=SA)

        spin = FindIC(primary,secondary)

        return spin


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
            self.parameters[parameter]=parameter_value
            setattr(self,parameter,parameter_value)
        if self.parameters['logQ']==True:self.convective_phase_lag=phase_lag(self.logQ)

        SecondaryAngmom=IntialSecondaryAngmom(self.interpolator,self.parameters)
        _logger.info('Seconary Initial Angular Momentum = {!r}'.format(repr(SecondaryAngmom())))

        binary_system=BinaryObjects(self.interpolator,self.parameters)

        primary=binary_system.create_star(self.primary_mass,dissipation=True)
        secondary=binary_system.create_star(self.secondary_mass,dissipation=True)
        binary=binary_system.create_binary_system(primary,secondary,secondary_angmom=SecondaryAngmom())

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
                        None,
                        timeout=3600)
            return binary

        

    def __init__(self,
                 interpolator,
                 parameters):

        self.interpolator=interpolator
        self.parameters=parameters

        for item,value in parameters.items():
            setattr(self,item,value)

        try:
            self.convective_phase_lag=phase_lag(self.logQ)
            _logger,info('Convective Phase lag = ',self.convective_phase_lag)
        except:
            self.parameter_logQ=True


