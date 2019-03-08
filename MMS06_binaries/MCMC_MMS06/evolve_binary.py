#!/usr/bin/env python3

"""Calculate orbital evolution of one of the MMS06 binaries."""

from argparse import ArgumentParser
import os.path
import pickle

import numpy
import astropy

from planetary_system_io import read_hatsouth_info
from reproduce_system import\
    find_evolution,\
    get_interpolator,\
    add_dissipation_cmdline,\
    get_poet_dissipation_from_cmdline
from stellar_evolution.manager import StellarEvolutionManager
from stellar_evolution import __path__ as poet_stellar_evolution_path

from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library

from fit_radial_velocities import fit_rv_data

def parse_command_line():
    """Parse the command line arguments into object with attributes."""

    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        'binary_info',
        help='The info file for the binary to evolve.'
    )
    add_dissipation_cmdline(parser)
    parser.add_argument(
        '--primary-stellar-evolution-interpolators',
        '--primary-interpolators',
        type=str,
        default=os.path.join(
            os.path.dirname(os.path.dirname(poet_stellar_evolution_path[0])),
            'stellar_evolution_interpolators'
        ),
        help='The path to a repository if stellar evolution interpolators.'
        ' Default: %(default)s'
    )
    parser.add_argument(
        '--secondary-stellar-evolution-interpolators',
        '--secondary-interpolators',
        type=str,
        default='stellar_evolution_interpolators',
        help='The path to a repository if stellar evolution interpolators.'
        ' Default: %(default)s'
    )
    parser.add_argument(
        '--secondary-track-path', '--secondary-tracks', '-t',
        type=str,
        default='MESA_like_tracks/binary6211',
        help='The directory containing the stellar evolution tracks to '
        'interpolate. Should have the same structure as the MESA tracks '
        'shipped with POET. Default: %(default)s'
    )
    parser.add_argument(
        '--eccentricity-expansion-coefficients', '--eccentricity-coefficients',
        type=str,
        default=os.path.join(
            os.path.dirname(os.path.dirname(poet_stellar_evolution_path[0])),
            'eccentricity_expansion_coef_O200.txt'
        ),
        help='The file with eccentricity expansion coefficients to use.'
    )
    parser.add_argument(
        '--initial-eccentricity', '--e0',
        type=float,
        default=0.0,
        help='The initial eccentricity to star the system with.'
    )
    parser.add_argument(
        '--initial-obliquity', '--obliquity',
        type=float,
        default=0.0,
        help='The initial obliquity between the orbit in which the secondary '
        'first forms and both zones of the primary star. Default: %(default)s'
    )
    parser.add_argument(
        '--disk-period', '--pdisk',
        type=float,
        default=10.0,
        help='The disk lock period to assume for the primary in days. '
        'Default: %(default)s.'
    )
    parser.add_argument(
        '--disk-lifetime',
        type=float,
        default=2e-3,
        help='The age when the secondary appears and the primary\'s spin is '
        'released in Gyr. Default: %(default)s.'
    )
    parser.add_argument(
        '--primary-angmom',
        type=float,
        nargs=3,
        default=(0.17, 2.78, 0.005),
        metavar=('<strength>',
                 '<saturation frequency>',
                 '<core-envelope copling timescale>'),
        help='The wind strength and saturation frequency, and the core-envelope'
        'coulping timescale defining the angular evolution for isolated stars '
        'to assume for the primary. Default: %(default)s.'
    )
    parser.add_argument(
        '--secondary-angmom',
        type=float,
        nargs=3,
        default=(0.01, 100, 0.005),
        metavar=('<strength>',
                 '<saturation frequency>',
                 '<core-envelope copling timescale>'),
        help='The wind strength and saturation frequency, and the core-envelope'
        'coulping timescale defining the angular evolution for isolated stars '
        'to assume for the secondary. Default: %(default)s.'
    )
    parser.add_argument(
        '--pickles',
        default='evolutions.pkl',
        help='A file containing pickles of past evolutions along with the '
        'configuration used to calculate them. Default: %(default)s.'
    )
    parser.add_argument(
        '--do-not-solve',
        action='store_true',
        default=False,
        help='If passed, no attempt is made to solve for the initial orbital '
        'period and/or eccentricity, instead the system parameters are assumed '
        'to be initial values.'
    )
    result = parser.parse_args()
    for component in ['primary', 'secondary']:
        for lgq_dependence in ['wtide', 'wspin']:
            num_args = len(
                getattr(result, '_'.join(('lgQ',
                                          component,
                                          lgq_dependence,
                                          'dependence')))
            )
            if num_args != 0 and (num_args < 3 or num_args % 2 == 0):
                parser.print_help()
                raise RuntimeError(
                    '--lgQ-ptide-dependence option requires a number of '
                    'argument that is at least 3 and odd!'
                )
    return result

def calculate_secondary_mass(primary_mass,
                             orbital_period,
                             rv_semi_amplitude,
                             eccentricity=0):
    """
    Calculate the mass of a secondary object for a given RV observations.

    All arguments must have the appropriate units set using the astropy module.

    Args:
        primary_mass:    The mass of the primary star.

        orbital_period:    The orbital period for the binary.

        rv_semi_amplitude:    The semi-amplitude of the radial velocity best-fit
            curve, complete with astropy units.

        eccentricity:    The eccentricity of the orbit.

    Returns:
        The mass of the secondary star, complete with astropy units.
    """

    mass_ratio_function = (
        (
            orbital_period * rv_semi_amplitude**3
            /
            #False positive
            #pylint: disable=no-member
            (2.0 * numpy.pi * astropy.constants.G)
            #pylint: enable=no-member
            *
            (1.0 - eccentricity**2)**1.5
        )
        /
        primary_mass
    ).to('')

    solutions = numpy.roots([1.0,
                             -mass_ratio_function,
                             -2.0 * mass_ratio_function,
                             -mass_ratio_function])
    mass_ratio = None
    for candidate_mass_ratio in solutions:
        if candidate_mass_ratio.imag == 0:
            assert mass_ratio is None
            mass_ratio = candidate_mass_ratio.real
    assert mass_ratio is not None
    return mass_ratio * primary_mass

def get_system_info(info_fname):
    """Return the system to evolve per the given info file."""

    info = read_hatsouth_info(info_fname)
    #False positive
    #pylint: disable=no-member
    orbital_parameters = fit_rv_data(info.RVdata)[0]
    info.rv_semi_amplitude = (orbital_parameters[0]
                              *
                              astropy.units.Unit('km/s'))
    info.eccentricity = orbital_parameters[1]
    info.Msecondary = calculate_secondary_mass(
        info.Mprimary,
        info.Porb,
        info.rv_semi_amplitude,
        info.eccentricity
    )
    #pylint: enable=no-member

    print(info.format())

    return info

def calculate_evolution(cmdline_args, info):
    """Calculate the evolution for a given command line configuration."""

    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        cmdline_args.eccentricity_expansion_coefficients.encode('ascii')
    )

    primary_interpolator = StellarEvolutionManager(
        cmdline_args.primary_stellar_evolution_interpolators
    ).get_interpolator_by_name('default')

    if (
            cmdline_args.secondary_stellar_evolution_interpolators
            ==
            cmdline_args.primary_stellar_evolution_interpolators
    ):
        secondary_interpolator = primary_interpolator
        #False positive
        #pylint: disable=no-member
        if (
                info.Msecondary.to(astropy.units.M_sun).value
                <
                secondary_interpolator.track_masses[0]
        ):
            info.Msecondary = (secondary_interpolator.track_masses[0]
                               *
                               astropy.units.M_sun)
        #pylint: enable=no-member
    else:
        secondary_interpolator = get_interpolator(
            cmdline_args.secondary_stellar_evolution_interpolators,
            cmdline_args.secondary_track_path
        )

    return find_evolution(
        info,
        (primary_interpolator, secondary_interpolator),
        get_poet_dissipation_from_cmdline(cmdline_args),
        initial_eccentricity=cmdline_args.initial_eccentricity,
        initial_inclination=cmdline_args.initial_obliquity,
        orbital_period_tolerance=1e-3,
        primary_wind_strength=cmdline_args.primary_angmom[0],
        primary_wind_saturation=cmdline_args.primary_angmom[1],
        primary_core_envelope_coupling_timescale=cmdline_args.primary_angmom[2],
        secondary_wind_strength=cmdline_args.secondary_angmom[0],
        secondary_wind_saturation=cmdline_args.secondary_angmom[1],
        secondary_core_envelope_coupling_timescale=(
            cmdline_args.secondary_angmom[2]
        ),
        #False positive
        #pylint: disable=no-member
        disk_period=(cmdline_args.disk_period * astropy.units.d),
        disk_dissipation_age=(cmdline_args.disk_lifetime
                              *
                              astropy.units.Gyr),
        #pylint: enable=no-member
        solve=(not cmdline_args.do_not_solve)
    )

def main(cmdline_args):
    """Calculate or use existing evolution per command line config."""

    system_info = get_system_info(cmdline_args.binary_info)

    evolution = None
    if os.path.exists(cmdline_args.pickles):
        with open(cmdline_args.pickles, 'rb') as pickle_file:
            while evolution is None:
                try:
                    pickled_cmdline = pickle.load(pickle_file)
                    if vars(pickled_cmdline) != vars(cmdline_args):
                        print(
                            '='*80
                            +
                            '\nCommand line mismatch:\n'
                            +
                            '\n'.join([
                                '%s: %s' % (repr(key), repr(value))
                                for key, value in vars(pickled_cmdline).items()
                            ])
                            +
                            '\nvs\n'
                            +
                            '\n'.join([
                                '%s: %s' % (repr(key), repr(value))
                                for key, value in vars(cmdline_args).items()
                            ])
                        )
                    pickled_info = pickle.load(pickle_file)
                    pickled_evolution = pickle.load(pickle_file)
                    if (
                            vars(pickled_cmdline) == vars(cmdline_args)
                            and
                            vars(pickled_info) == vars(system_info)
                    ):
                        print('Found existing evolution with specified '
                              'configuration')
                        evolution = pickled_evolution
                except EOFError:
                    break

    if evolution is None:
        evolution = calculate_evolution(cmdline_args, system_info)
        print('Evolution:\n' + evolution.format())
        with open(cmdline_args.pickles, 'ab') as pickle_file:
            pickle.dump(cmdline_args, pickle_file)
            pickle.dump(system_info, pickle_file)
            pickle.dump(evolution, pickle_file)

if __name__ == '__main__':
    main(parse_command_line())
