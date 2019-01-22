#!/usr/bin/env python3

"""Calculate orbital evolution of one of the MMS06 binaries."""

from argparse import ArgumentParser
from glob import glob
import os.path
import pickle

from matplotlib import pyplot
import numpy
import astropy

from planetary_system_io import read_hatsouth_info
from reproduce_system import find_evolution
from stellar_evolution.library_interface import MESAInterpolator
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
    parser.add_argument(
        '--lgQ-primary',
        type=float,
        default='6.0',
        help='The value of log10(Q*) to assume for the primary, at the '
        'reference tidal and spin periods if --lgQ-primary-wtide-dependence '
        'and/or --lgQ-primary-wspin-dependence is specified. '
        'Default: %(default)s.'
    )
    parser.add_argument(
        '--lgQ-secondary',
        type=float,
        default='6.0',
        help='The value of log10(Q*) to assume for the secondary, at the '
        'reference tidal and spin periods if --lgQ-secondary-wtide-dependence '
        'and/or --lgQ-secondary-wspin-dependence is specified. '
        'Default: %(default)s.'
    )
    parser.add_argument(
        '--lgQ-primary-wtide-dependence',
        nargs='+',
        type=float,
        default=[],
        metavar=('<powerlaw index> <break frequency> <powerlaw index>',
                 '<break frequency> <powerlaw index>'),
        help='Pass this argument to make lgQ of the primary depend on tidal '
        'period. At least three arguments must be passed: 1) the powerlaw index'
        ' for tidal frequencies below the first break, 2) the frequency '
        '[rad/day] where the first break occurs and 3) the powerlaw index after'
        ' the first break. Additional arguments must come in pairs, specifying '
        'more frequencies where breaks occur and the powerlaw indices for '
        'frequencies higher than the break.'
    )
    parser.add_argument(
        '--lgQ-secondary-wtide-dependence',
        nargs='+',
        type=float,
        default=[],
        metavar=('<powerlaw index> <break frequency> <powerlaw index>',
                 '<break frequency> <powerlaw index>'),
        help='Pass this argument to make lgQ of the secondary depend on tidal '
        'period. At least three arguments must be passed: 1) the powerlaw index'
        ' for tidal frequencies below the first break, 2) the frequency '
        '[rad/day] where the first break occurs and 3) the powerlaw index after'
        ' the first break. Additional arguments must come in pairs, specifying '
        'more frequencies where breaks occur and the powerlaw indices for '
        'frequencies higher than the break.'
    )
    parser.add_argument(
        '--lgQ-primary-wspin-dependence',
        nargs='+',
        type=float,
        default=[],
        metavar=('<powerlaw index> <break frequency> <powerlaw index>',
                 '<break frequency> <powerlaw index>'),
        help='Pass this argument to make lgQ of the primary depend on tidal '
        'period. At least three arguments must be passed: 1) the powerlaw index'
        ' for tidal frequencies below the first break, 2) the frequency '
        '[rad/day] where the first break occurs and 3) the powerlaw index after'
        ' the first break. Additional arguments must come in pairs, specifying '
        'more frequencies where breaks occur and the powerlaw indices for '
        'frequencies higher than the break.'
    )
    parser.add_argument(
        '--lgQ-secondary-wspin-dependence',
        nargs='+',
        type=float,
        default=[],
        metavar=('<powerlaw index> <break frequency> <powerlaw index>',
                 '<break frequency> <powerlaw index>'),
        help='Pass this argument to make lgQ of the secondary depend on tidal '
        'period. At least three arguments must be passed: 1) the powerlaw index'
        ' for tidal frequencies below the first break, 2) the frequency '
        '[rad/day] where the first break occurs and 3) the powerlaw index after'
        ' the first break. Additional arguments must come in pairs, specifying '
        'more frequencies where breaks occur and the powerlaw indices for '
        'frequencies higher than the break.'
    )
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
        '--pickles',
        default='evolutions.pkl',
        help='A file containing pickles of past evolutions along with the '
        'configuration used to calculate them. Default: %(default)s.'
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
                raise RuntimeError('--lgQ-ptide-dependence option requires a number of '
                                   'argument that is at least 3 and odd!')
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

def get_interpolator(stellar_evolution_interpolator_dir,
                     track_path):
    """Return the stellar evolution interpolator to use."""

    manager = StellarEvolutionManager(
        stellar_evolution_interpolator_dir
    )
    if not list(manager.get_suite_tracks()):
        manager.register_track_collection(
            track_fnames=glob(
                os.path.join(track_path, '*.csv')
            )
        )

    interpolator_args = dict(num_threads=1)
    interpolator_args['new_interp_name'] = 'custom'
    interpolator_args['nodes'] = {
        q: 0 for q in MESAInterpolator.quantity_list
    }
    interpolator_args['smoothing'] = {
        q: float('nan') for q in MESAInterpolator.quantity_list
    }
    return manager.get_interpolator(**interpolator_args)

def calculate_evolution(cmdline_config):
    """Calculate the evolution for a given command line configuration."""

    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        cmdline_config.eccentricity_expansion_coefficients.encode('ascii')
    )

    info = read_hatsouth_info(cmdline_config.binary_info)
    #False positive
    #pylint: disable=no-member
    info.rv_semi_amplitude = (fit_rv_data(info.RVdata)[0]
                              *
                              astropy.units.Unit('km/s'))
    info.Msecondary = calculate_secondary_mass(
        info.Mprimary,
        info.Porb,
        info.rv_semi_amplitude,
        info.eccentricity
    )
    #pylint: enable=no-member

    print(info.format())

    interpolator = (
        StellarEvolutionManager(
            cmdline_config.primary_stellar_evolution_interpolators
        ).get_interpolator_by_name('default'),
        get_interpolator(
            cmdline_config.secondary_stellar_evolution_interpolators,
            cmdline_config.secondary_track_path
        )
    )

    return find_evolution(
        info,
        interpolator,
        primary_lgq=cmdline_config.lgQ_primary,
        secondary_lgq=cmdline_config.lgQ_secondary,
        initial_eccentricity=cmdline_config.initial_eccentricity,
        orbital_period_tolerance=1e-3
    )

if __name__ == '__main__':
    cmdline_args = parse_command_line()


    evolution = None
    if os.path.exists(cmdline_args.pickles):
        with open(cmdline_args.pickles, 'rb') as pickle_file:
            while evolution is None:
                try:
                    if vars(pickle.load(pickle_file)) == vars(cmdline_args):
                        print('Found existing evolution with specified '
                              'configuration')
                        evolution = pickle.load(pickle_file)
                except EOFError:
                    break

    if evolution is None:
        evolution = calculate_evolution(cmdline_args)
        with open(cmdline_args.pickles, 'ab') as pickle_file:
            pickle.dump(cmdline_args, pickle_file)
            pickle.dump(evolution, pickle_file)

    print('Evolution: ' + evolution.format())

    pyplot.semilogx(
        #False positive
        #pylint: disable=no-member
        evolution.age,
        2.0 * numpy.pi / evolution.orbital_period,
        #pylint: enable=no-member
        label=r'$\Omega_{orb}$'
    )
    for component in ['primary', 'secondary']:
        pyplot.semilogx(
            #False positive
            #pylint: disable=no-member
            evolution.age,
            #pylint: enable=no-member
            (
                getattr(evolution, component + '_envelope_angmom')
                /
                getattr(evolution, component + '_iconv_star')
            ),
            label=component + r' $\Omega_{conv}$'
        )
        pyplot.semilogx(
            evolution.age,
            (
                getattr(evolution, component + '_core_angmom')
                /
                getattr(evolution, component + '_irad_star')
            ),
            label=component + r' $\Omega_{rad}$'
        )
    pyplot.twinx().semilogx(evolution.age,
                            evolution.eccentricity,
                            label='e')
    pyplot.legend()
    pyplot.show()
