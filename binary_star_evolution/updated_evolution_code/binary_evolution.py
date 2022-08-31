#!/usr/bin/env python3
import sys
from pathlib import Path
from directories import directories
import logging

home_dir=str(Path.home())
path=directories(home_dir)
sys.path.append(path.poet_path+'/PythonPackage')
sys.path.append(path.poet_path+'/scripts')

from orbital_evolution.transformations import phase_lag

from initial_secondary_angmom import IntialSecondaryAngmom
from nested_solver import InitialConditionSolver
from create_objects import BinaryObjects

_logger=logging.getLogger()

class Evolution:

    def calculate_intial_conditions(self):

        SecondaryAngmom=IntialSecondaryAngmom(self.interpolator,self.parameters)
        secondary_angmom=SecondaryAngmom()
        _logger.info('Seconary Initial Angular Momentum = {!r}'.format(repr(secondary_angmom)))

        binary_system=BinaryObjects(self.interpolator,self.parameters)

        primary=binary_system.create_star(self.primary_mass,dissipation=True)
        secondary=binary_system.create_star(self.secondary_mass,dissipation=True)


        FindIC=InitialConditionSolver(self.interpolator,
                                      self.parameters,
                                      secondary_angmom=secondary_angmom)

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
        # if 'logQ' in self.parameters['logQ']:self.convective_phase_lag=phase_lag(self.logQ)

        SecondaryAngmom=IntialSecondaryAngmom(self.interpolator,self.parameters)
        secondary_angmom=SecondaryAngmom()
        _logger.info('Seconary Initial Angular Momentum = {!r}'.format(secondary_angmom))

        binary_system=BinaryObjects(self.interpolator,self.parameters)

        primary=binary_system.create_star(self.primary_mass,dissipation=True)
        secondary=binary_system.create_star(self.secondary_mass,dissipation=True)
        binary=binary_system.create_binary_system(primary,
                                                  secondary,
                                                  initial_orbital_period=self.orbital_period,
                                                  initial_eccentricity=self.eccentricity,
                                                  secondary_angmom=secondary_angmom
                                                  )

        binary.evolve(self.age,
                      self.evolution_max_time_step,
                      self.evolution_precision,
                      None,
                      timeout=3600)



        final_state=binary.final_state()

        final_orbital_period=binary.orbital_period(final_state.semimajor)
        final_eccentricity=final_state.eccentricity
        print('final_eccntricity {} \nfinal orbital period {} \nfinal age {}'.format(final_eccentricity,final_orbital_period,final_state.age))

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
            _logger.info('Convective Phase lag = ',self.convective_phase_lag)
        except:
            self.parameter_logQ=True


