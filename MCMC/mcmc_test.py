#!/usr/bin/env python3

import matplotlib
matplotlib.use('TkAgg')

import sys
sys.path.append('../PythonPackage')
sys.path.append('../scripts')

import matplotlib.pyplot as plt
import scipy as sc
import numpy as np
from stellar_evolution.derived_stellar_quantities import TeffK
import scipy.interpolate
import scipy.linalg
import scipy.optimize
from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from orbital_evolution.binary import Binary
from orbital_evolution.transformations import phase_lag
from orbital_evolution.initial_condition_solver import InitialConditionSolver
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet
import numpy
from astropy import units, constants

pi=sc.pi


class struct:

    def __init__(self,**kwargs):

        for name,value in kwargs.items():
            setattr(self,name,value)


class MCMC:



    def mass_calculations(self):

        """returns primary and secondary star masses for the given parameters"""

    def create_planet(self ,mass = (constants.M_jup / constants.M_sun).to('')) :
        """Return a configured planet to use in the evolution."""

        planet = LockedPlanet(
            mass = mass,
            radius = (constants.R_jup / constants.R_sun).to('')
        )
        return planet


    def create_star(self, mass, wind=True):



        star = EvolvingStar(mass=mass,
                        metallicity=0.0,
                        wind_strength=0.17 if wind else 0.0,
                        wind_saturation_frequency=2.78,
                        diff_rot_coupling_timescale=5.0e-3,
                        interpolator=self.interpolator)

        star.set_dissipation(zone_index=0,
                         tidal_frequency_breaks=None,
                         spin_frequency_breaks=None,
                         tidal_frequency_powers=numpy.array([0.0]),
                         spin_frequency_powers=numpy.array([0.0]),
                         reference_phase_lag=self.convective_phase_lag)

        star.set_dissipation(zone_index=1,
                         tidal_frequency_breaks=None,
                         spin_frequency_breaks=None,
                         tidal_frequency_powers=numpy.array([0.0]),
                         spin_frequency_powers=numpy.array([0.0]),
                         reference_phase_lag=0.0)
        return star


    def create_binary_system(self,
                            primary,
                            secondary,
                            disk_dissipation_age,
                            secondary_angmom = None) :

        """Create a binary system to evolve from the given objects."""

        if isinstance(secondary, LockedPlanet) :
            secondary_config = dict(spin_angmom = numpy.array([0.0]),
                                inclination = None,
                                periapsis = None)
        else :
            secondary.select_interpolation_region(disk_dissipation_age)
            secondary_config = dict(spin_angmom = secondary_angmom,
                                inclination = numpy.array([0.0]),
                                periapsis = numpy.array([0.0]))

        secondary.configure(age = disk_dissipation_age,
                        companion_mass = primary.mass,
                        semimajor = self.initial_semimajor,
                        eccentricity = 0.0,
                        locked_surface = False,
                        zero_outer_inclination = True,
                        zero_outer_periapsis = True,
                        **secondary_config)

        if isinstance(secondary, EvolvingStar) :
            secondary.detect_stellar_wind_saturation()

        primary.select_interpolation_region(primary.core_formation_age())
        primary.detect_stellar_wind_saturation()


        binary = Binary(primary = primary,
                    secondary = secondary,
                    initial_orbital_period=self.porb_initial,
                    initial_semimajor = self.initial_semimajor,
                    initial_eccentricity = 0.0,
                    initial_inclination = 0.0,
                    disk_lock_frequency = self.disk_lock_frequency,
                    disk_dissipation_age = disk_dissipation_age,
                    secondary_formation_age = disk_dissipation_age)

        binary.configure(age = primary.core_formation_age(),
                     semimajor = float('nan'),
                     eccentricity = float('nan'),
                     spin_angmom = numpy.array([0.0]),
                     inclination = None,
                     periapsis = None,
                     evolution_mode = 'LOCKED_SURFACE_SPIN')

        return binary


    def binary_evolution(self):

        """returns current spin of the star given the parameters"""

        tdisk = 5e-3

        star = create_star(0.8, interpolator=interpolator, convective_phase_lag=0.0, wind=wind)
        planet = create_planet(1.0)

        binary = create_binary_system(star,
                                      planet,
                                      2.0 * numpy.pi / 3.0,
                                      10.0,
                                      tdisk)

        binary.evolve(tdisk, 1e-3, 1e-6, None)

        disk_state = binary.final_state()

        planet.delete()
        star.delete()
        binary.delete()

        primary = create_star(1.0, self.interpolator, self.convective_phase_lag, wind=wind)
        secondary = create_star(0.8, self.interpolator, self.convective_phase_lag, wind=wind)

        binary = create_binary_system(
            primary,
            secondary,
            2.0 * numpy.pi / 3.0,
            P_initial,
            tdisk,
            secondary_angmom=numpy.array([disk_state.envelope_angmom,
                                          disk_state.core_angmom])
        )

        binary.evolve(current_age, 1e-3, 1e-6, None)


        stellar_spin = 0.0

        primary.delete()
        secondary.delete()
        binary.delete()


        return  stellar_spin



    def acceptance_probability(self,current_values,proposed_values):

        prior_proposed =
        prior_current =


        likelihood_proposed =
        likelihood_current =

        p = (likelihood_proposed*prior_proposed)/(likelihood_current*prior_current)


        return p


    def values_proposed(self):

        for (name_obs,value_obs),(name_step,value_step) in zip(self.observation_data.items(),self.proposed_step.items()):

            current_value =

    def metropolis_hasting(self):

        """runs the metropolis hastrings algorithm for number of iterations given"""

        iteration_step =0






        current_values = initial_values

        f.open("accpted_parameters.txt")
        f.open("rejected_parameters.txt")

        # loop for number of iterations times

        while iteration_step < self.total_iterations:

            iteration_step = iteration_step + 1;



            # draw random value from the proposal funtion
            proposed_values = scipy.stats.norm.rvs(loc=current_values, scale=step_size)

            # calculate acceptance probablity
            p_acceptance = acceptance_probability(current_values, proposed_values)

            if p_acceptance > 1:

                current_values = proposed_values
                f.open("accepted_parameters.txt")
                f.write('step = %d\t', iteration_step, 'value = %d', current_values, '\n')
                f.close()

            else:

                f.open("rejected_parameters.txt")
                f.write('step = %d\t', iteration_step, 'value = %d', proposed_values, '\n')
                f.close()

                continue
        return none


    def __init__(self,interpolator,fixed_parameters,observation_data,logQ,proposed_step):


        """
        initialises parameters:
            -interpolator: stellar interpolator

            -binary_attributes:

                        -disk_dissipation_age

            -observables:

                        -Teff

                        -feh

                        -rvk

                        -inclination

            -total_iterations


        """

        self.interpolator  = interpolator
        self.fixed_parameters = fixed_parameters
        self.observation_data = observation_data
        self.convective_phase_lag = phase_lag(logQ)
        self.proposed_step = proposed_step

        self.update_parameters = dict()



if __name__ == '__main__':


    serialized_dir = "/home/kpenev/projects/git/poet/stellar_evolution_interpolators"
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

observation_data = dict(Teff = dict(vlaue = 1.0, sigma = 0.1),
                        feh = dict(value = 1.0, sigma = 0.1),
                        rvk = dict(value = 1.0, sigma = 0.1),
                        inclination = dict(value = 1.0, sigma = 0.1)
                    )

fixed_parameters = dict(disk_dissipation_age =  dict(value = 1.0, sigma = 0.1))

proposed_step = dict( Teff_step = 0.1,
                      feh_step = 0.1,
                      rvk_step = 0.1,
                      inclination = 0.1)



