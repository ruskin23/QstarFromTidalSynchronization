#!/usr/bin/env python3

import logging
import sys
import os
from pathlib import Path
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
import pickle

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

        initial_orbital_period=initial_conditions[0]
        initial_eccentricity=initial_conditions[1]



        _logger.info('\nTrying Porb_initial = {!r} , e_initial = {!r}'.format(initial_orbital_period,initial_eccentricity))

        
        if initial_eccentricity>0.80  or initial_eccentricity<0 or initial_orbital_period<0:
           return scipy.nan
        
        #if initial_eccentricity>0.80  or initial_eccentricity<0:
        #    _logger.warning('Encoutnered invalid value for eccentricity = {!r}'.format(initial_eccentricity))
        #    return -self.target_orbital_period,initial_eccentricity-self.target_eccentricity
        #if initial_orbital_period<0:
        #    _logger.warning('Encoutnered invalid value for orbtial period = {!r}'.format(initial_orbital_period))
        #    return initial_orbital_period-self.target_orbital_period,-self.target_eccentricity

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

        _logger.info(f'final_orbital_period = {self.final_orbital_period}\tfinal_eccentricity = {self.final_eccentricity}')


        if numpy.logical_or(numpy.isnan(self.final_orbital_period),numpy.isnan(self.final_eccentricity)):
            return scipy.nan
        #     evolution = binary.get_evolution()
        #     self.final_eccentricity,non_nan_index=check_last_nan(evolution.eccentricity)
        #     _logger.warning('Binary system was destroyed at age = {!r} Gyr'.format(evolution.age[non_nan_index]))
        #     self.delta_p=-self.target_orbital_period-self.target_age+evolution.age[non_nan_index]
        # else:
        #     self.delta_p=self.final_orbital_period-self.target_orbital_period
        # self.delta_e=self.final_eccentricity-self.target_eccentricity


        self.spin=(2*numpy.pi*binary.primary.envelope_inertia(final_state.age)/final_state.primary_envelope_angmom)

        binary.delete()

        _logger.info('Final Age = {!r}'.format(final_state.age))
        _logger.info('delta_p = {!r} , delta_e = {!r}'.format(self.delta_p,self.delta_e))
        _logger.info('Spin Period = %s',repr(self.spin))

        if self.function=='minimize':
            return numpy.sqrt(self.delta_p**2 + self.delta_e**2)
        else:
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

        initial_guess=[self.target_orbital_period,self.target_eccentricity]

        _logger.info('solving for p and e using function = {} and method {}'.format(self.function,self.method))

        Pguess=self.target_orbital_period
        e=self.target_eccentricity
        err_fun=numpy.nan
        n_1=0
        while numpy.isnan(err_fun):
            n_1=n_1+1
            Pguess=Pguess*n_1
            err_fun=self.initial_condition_errfunc([Pguess,e])
        simplex_1=(Pguess,e)

        err_fun=numpy.nan
        n_1=0
        while numpy.isnan(err_fun):
            n_1=n_1+1
            Pguess=2*Pguess*n_1
            err_fun=self.initial_condition_errfunc([Pguess,e])
        simplex_2=(Pguess,e)

        err_fun=numpy.nan
        
        while numpy.isnan(err_fun):
            e = e + (0.1 if e < 0.3 else -0.1)
            err_fun=self.initial_condition_errfunc([Pguess,e])
        simplex_3=(Pguess,e)

        initial_simplex=(simplex_1,simplex_2,simplex_3)
        _logger.info('Initial Simplex = {}'.format(repr(initial_simplex)))



        # try:
        #     if self.function=='least_squares':
        #         bounds=(numpy.array([0.0,0.0]), numpy.array([numpy.inf,0.8]))
        #         sol = optimize.least_squares(self.initial_condition_errfunc,
        #                         initial_guess,
        #                         bounds=bounds,
        #                         method=self.method
        #                         )
        #     if self.function=='root':
        #         sol = optimize.root(self.initial_condition_errfunc,
        #                         initial_guess,
        #                         method=self.method
        #                         )
        #     if self.function=='minimize':
        #         bounds=((0.0,numpy.inf), (0.0,0.8))
        #         sol = optimize.minimize(self.initial_condition_errfunc,
        #                         initial_guess,
        #                         bounds=bounds,
        #                         method=self.method
        #                         )


        #     initial_orbital_period_sol,initial_eccentricity_sol=sol.x
        # except:
        #     initial_orbital_period_sol,initial_eccentricity_sol=scipy.nan,scipy.nan
        #     self.spin=numpy.inf


        # _logger.info('Solver Results:')
        # _logger.info('Intial Orbital Period = {!r} , Initial Eccentricity = {!r}'.format(initial_orbital_period_sol,initial_eccentricity_sol))
        # _logger.info('Final Orbital Period = {!r} , Final Eccentricity = {!r}'.format(self.final_orbital_period,self.final_eccentricity))
        # _logger.info('Errors: delta_p = {!r} , delta_e = {!r}'.format(self.delta_p,self.delta_e))
        # _logger.info('Final Spin Period = {!r}'.format(self.spin))

        # results=dict()
        # results['p_initial']=initial_orbital_period_sol
        # results['e_initial']=initial_eccentricity_sol
        # results['p_final']=self.final_orbital_period
        # results['e_final']=self.final_eccentricity
        # results['delta_p']=self.delta_p
        # results['delta_e']=self.delta_e
        # results['spin']=self.spin


        # return results

