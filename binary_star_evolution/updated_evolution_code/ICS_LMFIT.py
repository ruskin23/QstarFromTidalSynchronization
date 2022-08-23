#!/usr/bin/env python3

from cmath import isnan
import logging
import sys
import os
from pathlib import Path
from unittest import result
from directories import directories

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
from basic_utils import Structure
from astropy import units, constants
from basic_utils import Structure
from create_objects import BinaryObjects

from math import pi
import scipy
from scipy import optimize
import numpy

from lmfit import minimize, Parameters, Parameter, report_fit

wsun = 0.24795522138
_logger = logging.getLogger()


def check_last_nan(a):
    for i in range(len(a)):
        val=a[len(a)-i-1]
        check=numpy.isnan(val)
        if check == False:
            return(val,len(a)-i-1)

class InitialConditionSolver:
    """Find initial orbital period and eccentricity which reproduce
        current orbital period and eccentricity of a given system """

    def initial_condition_errfunc(self,initial_conditions):
        """Error function which returns the difference between final values and intial values"""

        initial_orbital_period=initial_conditions['P_i'].value
        initial_eccentricity=initial_conditions['e_i'].value

        _logger.info('\nTrying Porb_initial = {!r} , e_initial = {!r}'.format(initial_orbital_period,initial_eccentricity))


        binary_system=BinaryObjects(self.interpolator,self.parameters)

        binary=binary_system.create_binary_system(self.primary,
                                                  self.secondary,
                                                  initial_orbital_period=initial_orbital_period,
                                                  initial_eccentricity=initial_eccentricity,
                                                  secondary_angmom=self.secondary_angmom)

        binary.evolve(
            self.age,
            self.evolution_max_time_step,
            self.evolution_precision,
            None,
            timeout=3600
        )

        final_state=binary.final_state()
        if final_state.age != self.target_age:
            _logger.warning('Evolution did not reach target age, crashed at age = {!r} Gyr'.format(final_state.age))
            assert(final_state.age==self.target_age)

        
        self.final_orbital_period=binary.orbital_period(final_state.semimajor)
        self.final_eccentricity=final_state.eccentricity

        if numpy.logical_or(numpy.isnan(self.final_orbital_period),numpy.isnan(self.final_eccentricity)):
            _logger.warning('Enountered NaN in final_orbital_period={!r} or final_eccentricity={!r}'.format(self.final_orbital_period,self.final_eccentricity))
            _logger.warning('Binary was destroyed')
            return numpy.array([scipy.nan,scipy.nan])
        #     evolution = binary.get_evolution()
        #     self.final_eccentricity,non_nan_index=check_last_nan(evolution.eccentricity)
        #     _logger.warning('Binary system was destroyed at age = {!r} Gyr'.format(evolution.age[non_nan_index]))
        #     self.delta_p=-self.target_orbital_period-self.target_age+evolution.age[non_nan_index]
        self.delta_p=self.final_orbital_period-self.target_orbital_period
        self.delta_e=self.final_eccentricity-self.target_eccentricity


        self.spin=(2*numpy.pi*binary.primary.envelope_inertia(final_state.age)/final_state.primary_envelope_angmom)

        binary.delete()

        _logger.info('final_orbital_period = {!r} , final_eccentricity = {!r}'.format(self.final_orbital_period,self.final_eccentricity))
        _logger.info('delta_p = {!r} , delta_e = {!r}'.format(self.delta_p,self.delta_e))
        _logger.info('Spin Period = %s',repr(self.spin))

        return numpy.array([self.delta_p,self.delta_e])



    def __init__(self,
                interpolator,
                parameters,
                evolution_max_time_step=1e-3,
                evolution_precision=1e-6,
                secondary_angmom=None):
        """
        Args:

            - interpolator

            - parameters
                a dictionary of all the parameters

            - evolution_max_time_step:
                The maximum timestep the evolution is allowed to make.

            - evolution_precision:
                The precision to require of the evolution.

            - secondary_angmom:
                The initial angular momentum of secondary star
        """

        self.interpolator=interpolator
        self.parameters=parameters
        for item,value in parameters.items():
            setattr(self,item,value)
        if 'logQ' in parameters: self.convective_phase_lag=phase_lag(self.logQ)
        else: self.convective_phase_lag=self.phase_lag_max


        self.evolution_max_time_step=evolution_max_time_step
        self.evolution_precision = evolution_precision
        self.secondary_angmom = secondary_angmom

        self.final_orbital_period,self.final_eccentricity=scipy.nan,scipy.nan
        self.delta_p,self.delta_e=scipy.nan,scipy.nan
        self.spin=scipy.nan
        self.initial_guess=None



    def calculate_good_initial_guess(self):
        
        P_guess=self.target_orbital_period
        e_guess=self.target_eccentricity
        params=Parameters()
        params.add('P_i',value=P_guess,min=0,max=100)
        params.add('e_i',value=e_guess,min=0,max=0.8)
        n=0
        while True:
            n=n+1
            dp,de=self.initial_condition_errfunc(params)
            if numpy.isnan(dp) or numpy.isnan(de):
                Pguess+=0.5
                continue
            else:
                self.initial_guess=[P_guess,e_guess]
                self.err_intial_guess=[dp,de]
                _logger.info('Found non Nan intial guess in {!r} try. dp={!r},de={!r}'.format(n,dp,de))
                break
            

    def __call__(self, primary, secondary):
        """
        Find initial conditions which reproduce the given system now.

        Args:
            - primary

            - secondary
        """

        self.primary = primary
        self.secondary = secondary

        self.target_age=self.age
        self.target_orbital_period=self.orbital_period
        self.target_eccentricity=self.eccentricity
        


        self.calculate_good_initial_guess()

        params=Parameters()
        params.add('P_i',value=self.initial_guess[0],min=0,max=60)
        params.add('e_i',value=self.initial_guess[1],min=0,max=0.8)

        method=self.method

        result=minimize(self.initial_condition_errfunc,params,method=method,nan_policy='omit')
        print(result)
        

        # _logger.info('Solver_Results:')
        # _logger.info('Intial_Orbital_Period={!r} , Initial_Eccentricity={!r}'.format(initial_orbital_period_sol,initial_eccentricity_sol))
        # _logger.info('Final_Orbital_Period={!r} , Final_Eccentricity={!r}'.format(self.final_orbital_period,self.final_eccentricity))
        # _logger.info('Errors: delta_p={!r} , delta_e={!r}'.format(self.delta_p,self.delta_e))
        _logger.info('Final_Spin_Period={!r}'.format(self.spin))

        return self.spin

