#!/usr/bin/env python3

"""Calculate orbital evolution of one of the MMS06 binaries."""

from argparse import ArgumentParser
from glob import glob
import os.path

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
        '--lgQ',
        type=float,
        default='6.0',
        help='The value of log10(Q*) to assume for both the primary and the '
        'secondary. Default: %(default)s.'
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
    return parser.parse_args()

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

if __name__ == '__main__':
    cmdline_args = parse_command_line()

    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        cmdline_args.eccentricity_expansion_coefficients.encode('ascii')
    )

    info = read_hatsouth_info(cmdline_args.binary_info)
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
            cmdline_args.primary_stellar_evolution_interpolators
        ).get_interpolator_by_name('default'),
        get_interpolator(
            cmdline_args.secondary_stellar_evolution_interpolators,
            cmdline_args.secondary_track_path
        )
    )

    evolution = find_evolution(info,
                               interpolator,
                               primary_lgq=cmdline_args.lgQ,
                               secondary_lgq=cmdline_args.lgQ,
                               initial_eccentricity=0.3)
