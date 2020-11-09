#!/usr/bin/env python3


import os
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
from multiprocessing import Pool

from math import pi
import scipy
from scipy import optimize
import numpy


wsun = 0.24795522138


class InitialConditionSolver:
    """Find initial orbital period and eccentricity which reproduce
        current orbital period and eccentricity of a given system """

    def initial_condition_errfunc(self,initial_p):
        """Error function which returns the difference between final values and intial values"""


        initial_orbital_period=initial_p
        #print('\nTrying Porb_initial = {} '.format(initial_orbital_period))

        binary_system=BinaryObjects(self.interpolator,self.parameters)

        binary=binary_system.create_binary_system(self.primary,
                                                  self.secondary,
                                                  initial_orbital_period=initial_orbital_period,
                                                  initial_eccentricity=self.initial_eccentricity,
                                                  secondary_angmom=self.secondary_angmom)

        binary.evolve(
            self.age,
            self.evolution_max_time_step,
            self.evolution_precision,
            None,
            timeout=3600
        )

        final_state=binary.final_state()
        assert(final_state.age==self.target_age)

        self.final_orbital_period=binary.orbital_period(final_state.semimajor)
        self.final_eccentricity=final_state.eccentricity

        self.delta_p=self.final_orbital_period-self.target_orbital_period
        self.delta_e=self.final_eccentricity-self.target_eccentricity

        if numpy.logical_or(numpy.isnan(self.delta_p),(numpy.isnan(self.delta_e))):
            print('Binary system was destroyed')
            raise ValueError

        self.spin=(2*numpy.pi*binary.primary.envelope_inertia(final_state.age)/final_state.primary_envelope_angmom)

        binary.delete()
        #print('delta_p = {} , delta_e = {}'.format(self.delta_p,self.delta_e))
        #print('Spin Period = ',self.spin)

        return self.delta_p

    def intial_guess(self):

        initial_p=self.target_orbital_period
        error_initial=self.initial_condition_errfunc(initial_p)

        if error_initial<0:
            p_min=initial_p
            while True:
                p_max=initial_p+1.0
                error_p=self.initial_condition_errfunc(p_max)
                if error_p>0:break
        else:
            p_max=initial_p
            while True:
                p_min=initial_p-1.0
                error_p=self.initial_condition_errfunc(p_min)
                if error_p<0:break

        return [p_min,p_max]

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

        self.evolution_max_time_step=evolution_max_time_step
        self.evolution_precision = evolution_precision
        self.secondary_angmom = secondary_angmom

        self.final_orbital_period=scipy.nan
        self.delta_p=scipy.nan
        self.spin=scipy.nan

    def solver(self,eccentricity):

        self.initial_eccentricity=eccentricity

        initial_guess=self.intial_guess()

        print('initial guess = ',initial_guess)

        print('solving for p')
        initial_orbital_period_sol = optimize.brentq(self.initial_condition_errfunc,
                                                    initial_guess[0],initial_guess[1],
                                                    xtol=1e-6
                                                    )

        return initial_orbital_period_sol,self.final_eccentricity,self.delta_p,self.delta_e,self.spin

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


        pool=Pool(processes=6)

        found_range=False
        found_solution=False
        delta_e_values=[]

        e_min=self.target_eccentricity
        e_max=0.43

        #use sleep function to see whats going on
        #use queue
        
        for k in range(2):            
            if k>0 and found_range==False:
                pool.terminate()
                pool.close()
                pool.join()
                break
            eccentricity_range=numpy.linspace(e_min,e_max,6)
            p=pool.imap(self.solver,eccentricity_range)

            for i,res in enumerate(p):

                porb_init,e_fin,delta_p,delta_e,spin=res
                porb_fin=self.target_orbital_period
                e_init=eccentricity_range[i]
                delta_e_values.append(delta_e)

                print('porb_init={} e_init={} porb_fin={} e_fin={}'.format(porb_init,e_init,porb_fin,e_fin))
                print('delta_p={} delta_e={}'.format(delta_p,delta_e))

                if i==0: 
                    if abs(delta_e)<1e-6:
                        print('ff found e solution = ',e_init)
                        found_solution=True
                        break
                    else:continue
                else:
                    if abs(delta_e_values[i])<1e-6:
                        print('found e solution = ',e_init)
                        found_solution=True
                        break

                    elif delta_e_values[i]*delta_e_values[i-1]<0:
                        e_min=eccentricity_range[i-1]
                        e_max=eccentricity_range[i]
                        print('found e_range = {} and {}'.format(e_min,e_max))
                        found_range=True
                        break
                    else:continue

        
        if found_solution==False:
            print('no solution found')

        print('Solver Results:')
        print('Intial Orbital Period = {} , Initial Eccentricity = {}'.format(porb_init,e_init))
        print('Final Orbital Period = {} , Final Eccentricity = {}'.format(porb_fin,e_fin))
        print('Errors: delta_p = {} , delta_e = {}'.format(delta_p,delta_e))
        print('Final Spin Period = {}'.format(spin))

