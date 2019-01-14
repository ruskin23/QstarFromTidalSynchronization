#!/usr/bin/env python3

import matplotlib

matplotlib.use('TkAgg')

# import sys
# sys.path.append('/Users/ruskinpatel/Desktop/Research/poet/PythonPackage')
# sys.path.append('/Users/ruskinpatel/Desktop/Research/poet/scripts')

import sys

sys.path.append('/home/kpenev/projects/git/poet/PythonPackage')
sys.path.append('/home/kpenev/projects/git/poet/scripts')

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


def create_star(mass, dissipation, interpolator, convective_phase_lag, wind=True):
    star = EvolvingStar(mass=mass,
                        metallicity=  -0.317892129759198,
                        wind_strength=0.17 if wind else 0.0,
                        wind_saturation_frequency=2.78,
                        diff_rot_coupling_timescale=5.0e-3,
                        interpolator=interpolator)

    if dissipation == 1:
        star.set_dissipation(zone_index=0,
                             tidal_frequency_breaks=None,
                             spin_frequency_breaks=None,
                             tidal_frequency_powers=numpy.array([0.0]),
                             spin_frequency_powers=numpy.array([0.0]),
                             reference_phase_lag=convective_phase_lag)

    return star


def create_binary_system(primary,
                         secondary,
                         disk_lock_frequency,
                         initial_orbital_period,
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

    if isinstance(secondary, EvolvingStar):
        secondary.detect_stellar_wind_saturation()

    primary.select_interpolation_region(primary.core_formation_age())
    primary.detect_stellar_wind_saturation()

    binary = Binary(primary=primary,
                    secondary=secondary,
                    initial_orbital_period=initial_orbital_period,
                    initial_eccentricity=0.0,
                    initial_inclination=0.0,
                    disk_lock_frequency=disk_lock_frequency,
                    disk_dissipation_age=disk_dissipation_age,
                    secondary_formation_age=disk_dissipation_age)

    secondary.configure(age=disk_dissipation_age,
                        companion_mass=primary.mass,
                        semimajor=binary.semimajor(initial_orbital_period),
                        eccentricity=0.0,
                        locked_surface=False,
                        zero_outer_inclination=True,
                        zero_outer_periapsis=True,
                        **secondary_config)

    binary.configure(age=primary.core_formation_age(),
                     semimajor=float('nan'),
                     eccentricity=float('nan'),
                     spin_angmom=numpy.array([0.0]),
                     inclination=None,
                     periapsis=None,
                     evolution_mode='LOCKED_SURFACE_SPIN')

    return binary


def plot_evolution(age,binary, wsat, style=dict(pcore='-b', penv='-g', score='m', senv='y')):
    """Calculate and plot the evolution of a properly constructed binary."""

    wsun = 0.24795522138  # 2*pi/25.34

    binary.evolve(age, 1e-3, 1e-6, None,True)

    evolution = binary.get_evolution()

    # print("wsun = ", wsun)
    # print("==   ", binary.secondary.core_inertia(evolution.age))

    wenv_secondary = (evolution.secondary_envelope_angmom / binary.secondary.envelope_inertia(evolution.age)) / wsun
    wcore_secondary = (evolution.secondary_core_angmom / binary.secondary.core_inertia(evolution.age)) / wsun
    wenv_primary = (evolution.primary_envelope_angmom / binary.primary.envelope_inertia(evolution.age)) / wsun
    wcore_primary = (evolution.primary_core_angmom / binary.primary.core_inertia(evolution.age)) / wsun

    quantities = ['primary_envelope', 'primary_core', 'secondary_envelope', 'secondary_core']

    orbitalfrequncy = binary.orbital_frequency(evolution.semimajor) / wsun

    pyplot.semilogx(evolution.age, wenv_primary, "-b", label='primary_envelope')
    pyplot.semilogx(evolution.age, wcore_primary, "-g", label='primary_core')
    pyplot.semilogx(evolution.age, wcore_secondary, "-m", label='secondary_core')
    pyplot.semilogx(evolution.age, wenv_secondary, "-r", label='secondary_envelope')

    pyplot.semilogx(evolution.age, orbitalfrequncy, "-k", label='orbital_frequency')
    pyplot.legend(loc='upper right')
    pyplot.ylim(top=30)
    pyplot.ylim(bottom=-20)
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
                                         'secondary_env_inertia'
                                         ]]))
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
    age = 2.6361887829110797
    primary_mass = 0.9735205066895405
    secondary_mass = 0.7522396868226867
    initial_disk_period =  1.4064967495370835 
    initial_orbital_period =5.2663825# 5.2663825
    


    star = create_star(secondary_mass, 1, interpolator=interpolator, convective_phase_lag=convective_phase_lag, wind=wind)
    planet = create_planet(1.0)

    binary = create_binary_system(star,
                                  planet,
                                  2.0 * numpy.pi / initial_disk_period,
                                  initial_orbital_period,
                                  tdisk)

    binary.evolve(tdisk, 1e-3, 1e-6, None)

    disk_state = binary.final_state()

    print("FINISHED PLANET-STAR EVOLUTION")

    planet.delete()
    star.delete()
    binary.delete()

    primary = create_star(primary_mass,1, interpolator, convective_phase_lag, wind=wind)
    secondary = create_star(secondary_mass, 1, interpolator, convective_phase_lag, wind=wind)
    # secondary = create_planet(1.0)

    print("Secondary_initial_angmom = ", numpy.array([disk_state.envelope_angmom, disk_state.core_angmom]))


    binary = create_binary_system(
        primary,
        secondary,
        2.0 * numpy.pi /initial_disk_period, 
        initial_orbital_period,
        tdisk,
        secondary_angmom=numpy.array([disk_state.envelope_angmom,
                                      disk_state.core_angmom])
    )

    print("BINARY STAR SYSTEM CREATED")

    # print(numpy.array([disk_state.envelope_angmom, disk_state.core_angmom]))

    evolution = plot_evolution(age,binary, wsat=2.78,
                               style=dict(orb='xr', core='xb', env='xg', sec_env=':c', sec_core=':m'))

    print("FINISHED BINARY STAR EVOLUTION")

    disk_state = binary.final_state()

    # print (disk_state.envelope_angmom)

    output_evolution(evolution, binary)

    print("ORBITAL_FREUQNCY = ", binary.orbital_frequency(evolution.semimajor))

    primary.delete()
    secondary.delete()
    binary.delete()


if __name__ == '__main__':
    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        b"eccentricity_expansion_coef.txt"
    )
    serialized_dir = '/home/kpenev/projects/git/poet/stellar_evolution_interpolators'

    # serialized_dir = '/Users/ruskinpatel/Desktop/Research/poet/stellar_evolution_interpolators'

    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    test_evolution(interpolator,  1.6799410609204806e-05, True)


