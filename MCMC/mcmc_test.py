
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


class data_set :

    def __init__(self,**kwargs):

        for name,value in kwargs.items():
            setattr(self,name,value)



def create_planet(mass = (constants.M_jup / constants.M_sun).to('')) :
    """Return a configured planet to use in the evolution."""

    planet = LockedPlanet(
        mass = mass,
        radius = (constants.R_jup / constants.R_sun).to('')
    )
    return planet


def create_star(mass, interpolator, convective_phase_lag, wind=True):
    star = EvolvingStar(mass=mass,
                        metallicity=0.0,
                        wind_strength=0.17 if wind else 0.0,
                        wind_saturation_frequency=2.78,
                        diff_rot_coupling_timescale=5.0e-3,
                        interpolator=interpolator)

    star.set_dissipation(zone_index=0,
                         tidal_frequency_breaks=None,
                         spin_frequency_breaks=None,
                         tidal_frequency_powers=numpy.array([0.0]),
                         spin_frequency_powers=numpy.array([0.0]),
                         reference_phase_lag=convective_phase_lag)

    star.set_dissipation(zone_index=1,
                         tidal_frequency_breaks=None,
                         spin_frequency_breaks=None,
                         tidal_frequency_powers=numpy.array([0.0]),
                         spin_frequency_powers=numpy.array([0.0]),
                         reference_phase_lag=0.0)
    return star

def create_binary_system(primary,
                         secondary,
                         disk_lock_frequency,
                         initial_semimajor,
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
                        semimajor = initial_semimajor,
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
                    initial_semimajor = initial_semimajor,
                    initial_eccentricity = 0.0,
                    initial_inclination = 0.0,
                    disk_lock_frequency = disk_lock_frequency,
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



def initial_condition(current_age,current_Porb,current_Psurf,primary,secondary):

    find_ic = InitialConditionSolver(disk_dissipation_age = 5e-3,
                                     evolution_max_time_step = 1e-2)

    target = Structure(age = current_age,
                       Porb = current_Porb,
                       Psurf = current_Psurf,
                       planet_formation_age = 5e-3)


    initial_porb, initial_psurf = find_ic(target = target,
                                          star = primary,
                                          planet = secondary)


    return  initial_porb


def binary_evolution(current_age,current_Porb, current_Psurf,interpolator, convective_phase_lag, wind):

    """run evolution for binary system """

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

    primary = create_star(1.0, interpolator, convective_phase_lag, wind=wind)
    secondary = create_star(0.8, interpolator, convective_phase_lag, wind=wind)

    P_initial =  initial_condition(current_age, current_Porb, current_Psurf, primary, secondary)

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

    disk_state = binary.final_state()

    primary.delete()
    secondary.delete()
    binary.delete()

    return  disk_state.envelope_angmom



def acceptance_probability(values_current,values_proposed):

    p = (likelihood(values_proposed,sigma)*prior.sample_value(values_proposed))/(likelihood(values_current,sigma)*prior.sample_value(values_current))

    return p


def metropolis_hasting(initial_value,step_size,total_iterations):

    iteration_step =0;

    value_current = initial_value


    #loop for number of iterations times

    while True:

        iteration_step = iteration_step + 1;

        if iteration_step > total_iterations:

            break

        #draw random value from the proposal funtion
        value_proposed = scipy.stats.norm.rvs(loc = value_current,scale = step_size)


        #calculate acceptance probablity
        p_acceptance = acceptance_probability(value_current,value_proposed)

        if p_acceptance > 1 :

            value_current = value_proposed
            f.open("accpted_parameters.txt")
            f.write('step = %d\t',iteration_step, 'value = %d',value_current,'\n')
            f.close()

        else:

            f.open("rejected_parameters.txt")
            f.write('step = %d\t',iteration_step, 'value = %d',value_current,'\n')
            f.close()

            continue







