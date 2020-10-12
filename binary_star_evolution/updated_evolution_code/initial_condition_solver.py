#!/usr/bin/env python3

import shelve

import sys

sys.path.append('/home/ruskin/projects/poet/PythonPackage')
sys.path.append('/home/ruskin/projects/poet/scripts')


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

wsun = 0.24795522138


class InitialConditionSolver:
    """Find initial conditions which reproduce a given system now."""

    def initial_condition_errfunc(self,initial_conditions): 
        """Error function which returns the difference between final values and intial values"""

        initial_orbital_period=initial_conditions[0]
        initial_eccentricity=initial_conditions[1]

        print('\nTrying Porb_initial = {} , e_initial = {}'.format(initial_orbital_period,initial_eccentricity))

        if initial_eccentricity>0.45 or initial_orbital_period<0 or initial_eccentricity<0:
            print('Invalid values')
            return scipy.nan,scipy.nan
       
        binary_system=BinaryObjects(self.interpolator,self.parameters)
        
        binary=binary_system.create_binary_system(self.primary,self.secondary,secondary_angmom=self.secondary_angmom)
        
        binary.evolve(
            self.age, 
            self.evolution_max_time_step, 
            self.evolution_precision, 
            None,
            timeout=600
        )
                        
        final_state=binary.final_state()
        assert(final_state.age==self.age) 

        self.final_orbital_period=binary.orbital_period(final_state.semimajor)
        self.final_eccentricity=final_state.eccentricity

        self.delta_p=self.final_orbital_period-initial_orbital_period
        self.delta_e=self.final_eccentricity-initial_eccentricity

        if numpy.logical_or(numpy.isnan(self.delta_p),(numpy.isnan(self.delta_e))):
            print('Binary system was destroyed')
            raise ValueError

        self.spin=(2*numpy.pi*binary.primary.envelope_inertia(final_state.age)/final_state.primary_envelope_angmom)

        binary.delete()

        print('delta_p = {} , delta_e = {}'.format(self.delta_p,self.delta_e))
        print('Spin Period = ',self.spin)

        return self.delta_p,self.delta_e
   
    

    def __init__(self,
                interpolator,
                parameters,
                evolution_max_time_step=1e-3,
                evolution_precision=1e-6,
                secondary_angmom=None):
        """
        Initialize the intiial conditino solver object

        Args:

            - interpolator
            
            - parameters

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

        self.evolution_max_time_step=evolution_max_time_step
        self.evolution_precision = evolution_precision
        self.secondary_angmom = secondary_angmom

    def __call__(self, primary, secondary):
        """
        Find initial conditions which reproduce the given system now.

        Args:
            - primary

            - secondary
        """

        self.primary = primary
        self.secondary = secondary

        initial_guess=[self.orbital_period,self.eccentricity]


        print('solving for p and e')
        sol = optimize.root(self.initial_condition_errfunc,
                            initial_guess,
                            method='lm',
                            tol=1e-6,
                            options={'xtol':1e-6,
                                     'ftol':1e-6}
                            )

        initial_orbital_period_sol,initial_eccentricity_sol=sol.x

        print('Solver Results:')
        print('Intial Orbital Period = {} , Initial Eccentricity = {}'.format(initial_orbital_period_sol,initial_eccentricity_sol))
        print('Final Orbital Period = {} , Final Eccentricity = {}'.format(self.final_orbital_period,self.final_eccentricity))
        print('Errors: delta_p = {} , delta_e = {}'.format(self.delta_p,self.delta_e))
        print('Final Spin Period = {}'.format(self.spin))

