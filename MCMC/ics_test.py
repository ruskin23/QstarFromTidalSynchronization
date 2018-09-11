#!/usr/bin/env python3

import matplotlib
matplotlib.use('TkAgg')

import sys
sys.path.append('../PythonPackage')
sys.path.append('../scripts')

from matplotlib import pyplot
from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from orbital_evolution.binary import Binary
from orbital_evolution.transformations import phase_lag
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.initial_condition_solver import InitialConditionSolver
import numpy
from astropy import units, constants

wsun = 2.0 * numpy.pi / 25.34


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


def test_ic_solver(interpolator) :
    """Find initial condition to reproduce some current state and plot."""

    find_ic = InitialConditionSolver(disk_dissipation_age = 5e-3,
                                     evolution_max_time_step = 1e-2)
    target = Structure(age = 5.0,
                       Porb = 3.0,
                       Wdisk = 5.0,
                       planet_formation_age = 5e-3)

    primary = create_star(1.0, interpolator, convective_phase_lag, wind = wind)
    secondary = create_star(0.8, interpolator, convective_phase_lag, wind = wind)

    initial_porb, initial_psurf = find_ic(target = target,
                                          star = primary,
                                          planet = secondary)

    print('IC: Porb0 = %s, P*0 = %s' % (repr(initial_porb),
                                        repr(initial_psurf)))


if __name__ == '__main__' :

    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        b"eccentricity_expansion_coef.txt"
    )

    serialized_dir = '/home/kpenev/projects/git/poet/stellar_evolution_interpolators'
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    test_ic_solver(interpolator)
