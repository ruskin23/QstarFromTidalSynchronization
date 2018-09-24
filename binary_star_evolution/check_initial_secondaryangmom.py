#!/usr/bin/env python3

#checking the effects of value of initial semimajor axis on calculations of initial angular momentum of secondary
#RESULT = initial secondary angular momentum is unaffected by the choice of initial semi major axis, as the evolution stops
#very early because of the small value of disk dissipation age

import matplotlib

matplotlib.use('TkAgg')

import sys

sys.path.append('/Users/ruskinpatel/Desktop/Research/poet/PythonPackage')
sys.path.append('/Users/ruskinpatel/Desktop/Research/poet/scripts')

from matplotlib import pyplot
from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as \
    orbital_evolution_library
from orbital_evolution.binary import Binary
from orbital_evolution.transformations import phase_lag
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet
import numpy
from astropy import units, constants

wsun = 0.24795522138


def create_planet(mass=(constants.M_jup / constants.M_sun).to('')):
    """Return a configured planet to use in the evolution."""

    planet = LockedPlanet(
        mass=mass,
        radius=(constants.R_jup / constants.R_sun).to('')
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
                         secondary_angmom=None):
    """Create a binary system to evolve from the given objects."""

    if isinstance(secondary, LockedPlanet):
        secondary_config = dict(spin_angmom=numpy.array([0.0]),
                                inclination=None,
                                periapsis=None)
    else:
        secondary.select_interpolation_region(disk_dissipation_age)
        secondary_config = dict(spin_angmom=secondary_angmom,
                                inclination=numpy.array([0.0]),
                                periapsis=numpy.array([0.0]))

    secondary.configure(age=disk_dissipation_age,
                        companion_mass=primary.mass,
                        semimajor=initial_semimajor,
                        eccentricity=0.0,
                        locked_surface=False,
                        zero_outer_inclination=True,
                        zero_outer_periapsis=True,
                        **secondary_config)

    if isinstance(secondary, EvolvingStar):
        secondary.detect_stellar_wind_saturation()

    primary.select_interpolation_region(primary.core_formation_age())
    primary.detect_stellar_wind_saturation()

    binary = Binary(primary=primary,
                    secondary=secondary,
                    initial_semimajor=initial_semimajor,
                    initial_eccentricity=0.0,
                    initial_inclination=0.0,
                    disk_lock_frequency=disk_lock_frequency,
                    disk_dissipation_age=disk_dissipation_age,
                    secondary_formation_age=disk_dissipation_age)

    binary.configure(age=primary.core_formation_age(),
                     semimajor=float('nan'),
                     eccentricity=float('nan'),
                     spin_angmom=numpy.array([0.0]),
                     inclination=None,
                     periapsis=None,
                     evolution_mode='LOCKED_SURFACE_SPIN')

    return binary


def plot_evolution(binary, wsat, style=dict(core='-b', env='-g')):
    """Calculate and plot the evolution of a properly constructed binary."""

    wsun = 0.24795522138  # 2*pi/25.34

    binary.evolve(10.0, 1e-3, 1e-6, None)

    evolution = binary.get_evolution()

    print("wsun = ", wsun)
    # print("==   ", binary.secondary.core_inertia(evolution.age))

    wenv = (evolution.secondary_envelope_angmom / binary.secondary.envelope_inertia(evolution.age)) / wsun
    wcore = (evolution.secondary_core_angmom / binary.secondary.core_inertia(evolution.age)) / wsun

    pyplot.semilogx(evolution.age, wenv, style['env'])
    pyplot.semilogx(evolution.age, wcore, style['core'])

    pyplot.semilogx(evolution.age, binary.orbital_frequency(evolution.semimajor), "-k")

    pyplot.show()

    return evolution


def output_evolution(evolution, binary):
    """Write the given evolution to stdout organized in columns."""

    quantities = list(
        filter(lambda q: q[0] != '_' and q != 'format', dir(evolution))
    )
    print(' '.join(['%30s' % q for q in quantities]), end=' ')
    print(' '.join(['%30s' % q for q in ['primary_core_inertia',
                                         'primary_env_inertia',
                                         'secondary_core_inertia',
                                         'secondary_env_inertia']]))
    for i in range(len(evolution.age)):
        print(' '.join(['%30s' % repr(getattr(evolution, q)[i])
                        for q in quantities]), end=' ')
        age = evolution.age[i]
        print(' '.join([
            '%30.16e' % q for q in [
                binary.primary.core_inertia(age),
                binary.primary.envelope_inertia(age),
                binary.secondary.core_inertia(age),
                binary.secondary.envelope_inertia(age)
            ]
        ]))


def test_evolution(interpolator, convective_phase_lag, wind):
    """run evolution for binary system """

    tdisk = 5e-3

    star = create_star(0.8, interpolator=interpolator, convective_phase_lag=0.0, wind=wind)
    planet = create_planet(1.0)

    initial_semimajor = numpy.linspace(1.0,20.0,100)

    for a in initial_semimajor:

        binary = create_binary_system(star,
                                     planet,
                                      2.0 * numpy.pi / 3.0,
                                      10.0,
                                      tdisk)

        binary.evolve(tdisk, 1e-3, 1e-6, None)

        disk_state = binary.final_state()

        secondary_angmom = numpy.array([disk_state.envelope_angmom,
                                    disk_state.core_angmom])

        binary.delete()

        print (secondary_angmom)




    planet.delete()
    star.delete()


    print (secondary_angmom)

if __name__ == '__main__':
    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        b"eccentricity_expansion_coef.txt"
    )
    # serialized_dir = '/home/kpenev/projects/git/poet/stellar_evolution_interpolators'

    serialized_dir = '/Users/ruskinpatel/Desktop/Research/poet/stellar_evolution_interpolators'

    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    test_evolution(interpolator, phase_lag(6.0), True)


