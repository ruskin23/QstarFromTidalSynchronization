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

wsun = 0.24795522138
pi=sc.pi


class MCMC:

    def create_planet(self,mass = (constants.M_jup / constants.M_sun).to('')) :
        """Return a configured planet to use in the evolution."""

        planet = LockedPlanet(
            mass = mass,
            radius = (constants.R_jup / constants.R_sun).to('')
        )
        return planet


    def create_star(self ,mass, wind=True):

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


    def initial_condition(self):

        find_ic = InitialConditionSolver(disk_dissipation_age = 5e-3,
                                        evolution_max_time_step = 1e-2)
        target = Structure(age = self.age,
                           Porb = self.Porb,
                           Psurf = Psurf,
                           planet_formation_age = 5e-3)

        primary = create_star(1.0, self.interpolator, self.convective_phase_lag, wind=wind)
        secondary = create_star(0.8, self.interpolator, self.convective_phase_lag, wind=wind)

        initial_porb, initial_psurf = find_ic(target = target,
                                          star = primary,
                                          planet = secondary)


        return  initial_porb


def acceptance_probability(values_current,values_proposed):

    dist =  normal_distribution(mu,sigma)

    p = (L(values_proposed)*dist.sample_value(values_proposed))/(L(values_current)*dist.sample_value(values_current))

    return p


def metropolis_hasting(initial_value,step_size,total_iterations):

    iteration_step =0 ;

    value_current = initial_value

    f.open("accpted_parameters.txt")
    f.open("rejected_parameters.txt")

    #loop for number of iterations times

    while True:

        iteration_step = iteration_step + 1;

        #draw random value from the proposal funtion
        value_proposed = scipy.stats.norm.rvs(loc = value_current,scale = step_size)


        #calculate acceptance probablity
        p_acceptance = acceptance_probability(value_current,value_proposed)

        if p_acceptance > 1 :

            value_current = value_proposed


        else:

            continue

        if iteration_step < total_iterations:

            break





