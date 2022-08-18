#!/usr/bin/env python3

from cmath import isnan
import logging
import sys
import os
from pathlib import Path
from turtle import numinput
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

        if self.initial_guess is not None:
            if initial_orbital_period==self.initial_guess[0] and initial_eccentricity==self.initial_guess[1]:
                return numpy.array(self.err_intial_guess)

        if numpy.isnan(initial_orbital_period) or numpy.isnan(initial_eccentricity):
            _logger.warning('Solver using NaN as initial values.')
            raise ValueError('Solution cant be found as solver is trying NaN as initial values')

        if initial_eccentricity>0.80  or initial_eccentricity<0 or initial_orbital_period<0:
            _logger.warning('Encoutnered invalid initial values, returning NaN')
            if self.function=='minimize': return scipy.nan
            else: return numpy.array([scipy.nan,scipy.nan])

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
            if self.function=='minimize': return scipy.nan
            else: return numpy.array([scipy.nan,scipy.nan])
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

        if self.function=='minimize':
            err_fun =  numpy.sqrt(self.delta_p**2 + self.delta_e**2)
            _logger.info('Error Function = {!r}'.format(err_fun))
            return err_fun
        else: return numpy.array([self.delta_p,self.delta_e])



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


    def calculate_initial_simplex(self):

        Pguess=self.target_orbital_period
        e=self.target_eccentricity
        _logger.info('\nCalculating first simplex using Pguess = {} e =  {}'.format(repr(Pguess),repr(e)))
        err_fun=self.initial_condition_errfunc([Pguess,e])
        if numpy.isnan(err_fun):
            _logger.warning('Encontered nan values.')
            n=1
            while numpy.isnan(err_fun):
                n=n+1
                Pguess=Pguess*n
                err_fun=self.initial_condition_errfunc([Pguess,e])
        simplex_1=[Pguess,e]
        _logger.info('Found first simplex values at Pguess = {} e = {}'.format(repr(Pguess),repr(e)))

        Pguess=2*Pguess
        _logger.info('\nCalculating second simplex using Pguess = {} e =  {}'.format(repr(Pguess),repr(e)))
        err_fun=self.initial_condition_errfunc([Pguess,e])
        if numpy.isnan(err_fun):
            _logger.warning('Encontered nan values.')
            n=1
            while numpy.isnan(err_fun):
                n=n+1
                Pguess=Pguess*n
                err_fun=self.initial_condition_errfunc([Pguess,e])
        simplex_2=[Pguess,e]
        _logger.info('Found second simplex values at Pguess = {} e = {}'.format(repr(Pguess),repr(e)))


        e = e + 0.1 if e < 0.3 else -0.1
        _logger.info('\nCalculating third simplex using Pguess = {} e =  {}'.format(repr(Pguess),repr(e)))
        err_fun=self.initial_condition_errfunc([Pguess,e])
        if numpy.isnan(err_fun):
            _logger.warning('Encontered nan values.')
            while numpy.isnan(err_fun):
                e = e + (0.1 if e < 0.3 else -0.1)
                err_fun=self.initial_condition_errfunc([Pguess,e])
        simplex_3=[Pguess,e]
        _logger.info('Found second simplex values at Pguess = {} e = {}'.format(repr(Pguess),repr(e)))


        initial_simplex=numpy.array([simplex_1,simplex_2,simplex_3])
        _logger.info('\nInitial Simplex = {}'.format(repr(initial_simplex)))

        return initial_simplex

    def calculate_good_initial_guess(self):

        Pguess=self.target_orbital_period
        e=self.target_eccentricity
        while True:
            dp,de=self.initial_condition_errfunc([Pguess,e])
            if numpy.isnan(dp) or numpy.isnan(de):
                Pguess+=0.5
                continue
            else:
                self.initial_guess=[Pguess,e]
                self.err_intial_guess=[dp,de]
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
        

        _logger.info('\nSolving for p and e using function = {} and method {}'.format(self.function,self.method))


        if self.function=='minimize':
            self.initial_guess=[self.target_orbital_period,self.target_eccentricity]
            initial_simplex=self.calculate_initial_simplex()
        else:
            self.calculate_good_initial_guess()
        

        try:
            if self.function=='least_squares':
                bounds=(numpy.array([0.0,0.0]), numpy.array([numpy.inf,0.8]))
                sol = optimize.least_squares(self.initial_condition_errfunc,
                                self.initial_guess,
                                bounds=bounds,
                                method=self.method
                                )
            if self.function=='root':
                sol = optimize.root(self.initial_condition_errfunc,
                                self.initial_guess,
                                method=self.method,
                                options={'xtol': 0.0, 'ftol': 1e-05}
                                )
            if self.function=='minimize':
                bounds=((0.0,numpy.inf), (0.0,0.8))
                sol = optimize.minimize(self.initial_condition_errfunc,
                                self.initial_guess,
                                method=self.method,
                                bounds=bounds,
                                options={'initial_simplex' : initial_simplex}
                                )


            initial_orbital_period_sol,initial_eccentricity_sol=sol.x
        except:
            initial_orbital_period_sol,initial_eccentricity_sol=scipy.nan,scipy.nan
            self.spin=numpy.inf


        _logger.info('Solver_Results:')
        _logger.info('Intial_Orbital_Period={!r} , Initial_Eccentricity={!r}'.format(initial_orbital_period_sol,initial_eccentricity_sol))
        _logger.info('Final_Orbital_Period={!r} , Final_Eccentricity={!r}'.format(self.final_orbital_period,self.final_eccentricity))
        _logger.info('Errors: delta_p={!r} , delta_e={!r}'.format(self.delta_p,self.delta_e))
        _logger.info('Final_Spin_Period={!r}'.format(self.spin))

        if numpy.isfinite(self.spin):
            sum_of_squares=numpy.sqrt(self.delta_p**2+self.delta_e**2)
            if sum_of_squares>1e-3:
                _logger.info('Error Margins not good enough. Setting Spin=inf')
                self.spin=numpy.inf
            if numpy.isnan(self.spin):
                _logger.info('Spin=NaN. Setting Spin=inf')
                self.spin=numpy.inf
        if numpy.isnan(self.spin):
            _logger.warning('Spin=Nan after solver. Setting Spin=inf')
            self.spin=numpy.inf

        return self.spin

