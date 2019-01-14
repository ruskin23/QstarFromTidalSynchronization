#!/usr/bin/env python3
"""Fit or the amplitude, offset and shift describing an RV orbit."""

from argparse import ArgumentParser

import scipy
import scipy.optimize

def parse_command_line():
    """Parse the command line to object with attributes."""

    parser = ArgumentParser(
        description='Fit for the radial velocity semi-amplitude from a set of '
        'RV points. Assumes circular orbit and that the x-axis is the orbital '
        'phase.'
    )
    parser.add_argument(
        'rv_fname',
        help='The name of the file containing the RV measurements to fit.'
    )
    return parser.parse_args()

def residuals(params, rv_data):
    """
    Return the residuals of the RV estimate with the given parameters.

    Args:
        params:    The amplitude, phase and offset of the RV dependence on
            phase.

        rv_data:    The available radial velocity measurements.

    Returns:
        scipy array:
            The residuals of the RV function per the given parameters and RV
            data.
    """

    amplitude, phase, offset = params
    return (
        amplitude * scipy.cos(2.0 * scipy.pi * rv_data[:, 0] + phase)
        +
        offset
        -
        rv_data[:, 1]
    )

def fit_rv_data(rv_data_fname):
    """Fit or the amplitude, offset and shift describing an RV orbit."""

    rv_data = scipy.genfromtxt(rv_data_fname)
    (amplitude, phase, offset), flag = scipy.optimize.leastsq(residuals,
                                                              [1.0, 0.0, 0.0],
                                                              args=(rv_data,))
    assert(flag in [1, 2, 3, 4])
    return amplitude, phase, offset

if __name__ == '__main__':
    cmdline_args = parse_command_line()
    fit_results = fit_rv_data(cmdline_args.rv_fname)
    print('Amplitude: ' + repr(fit_results[0]))
