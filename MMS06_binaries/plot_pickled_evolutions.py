#!/usr/bin/env python3

"""Plot the evolutions contained in a pickle file created by evolve_binary."""

import pickle
from argparse import ArgumentParser

from matplotlib import pyplot
import numpy

def parse_command_line():
    """Parse the command line arguments into object with attributes."""

    parser = ArgumentParser(description=__doc__)

    parser.add_argument(
        '--pickle',
        default='evolutions.pkl',
        help='The pickle file to open. Default: %(default)s.'
    )
    parser.add_argument(
        '--to-file',
        nargs=2,
        action='append',
        metavar=('<evolution index>', '<filename>'),
        default=[],
        help='If this argument is specified, only the selected evolutions are '
        ',plotted and the resulting plots are saved to files with the '
        'specified names rather than displayed. Can be used multiple times.'
    )
    parser.add_argument(
        '--exclude-core-spin',
        action='store_true',
        help='If passed, the plots will not include the evolution of the spin '
        'of the radiative core(s) of the components.'
    )
    return parser.parse_args()

def plot_evolution(evolution, cmdline_config=None, exclude_core_spin=False):
    """Plot the evolution and show it."""

    w_plot_include = evolution.age > 0.002
    w_axis = pyplot.gca()
    e_axis = pyplot.twinx()

    plot_ages = evolution.age[w_plot_include] * 1e3
    if cmdline_config is not None:
        plot_ages -= getattr(cmdline_config, 'disk_lifetime', 2e-3) * 1e3

    w_axis.semilogx(
        #False positive
        #pylint: disable=no-member
        plot_ages,
        2.0 * numpy.pi / evolution.orbital_period[w_plot_include],
        #pylint: enable=no-member
        label=r'$\Omega_{orb}$',
        linewidth=5
    )
    for component in ['primary', 'secondary']:
        component_label = 1 if component == 'primary' else 2
        w_axis.semilogx(
            plot_ages,
            (
                getattr(evolution, component + '_envelope_angmom')
                /
                getattr(evolution, component + '_iconv')
            )[w_plot_include],
            label=(r' $\Omega_%s^{conv}$' % component_label),
            linewidth=3
        )
        if not exclude_core_spin:
            w_axis.semilogx(
                plot_ages,
                (
                    getattr(evolution, component + '_core_angmom')
                    /
                    getattr(evolution, component + '_irad')
                )[w_plot_include],
                label=(r' $\Omega_%s^{rad}$' % component_label),
                linewidth=3
            )

    e_axis.semilogx(plot_ages,
                    evolution.eccentricity[w_plot_include],
                    '-k',
                    label='e',
                    linewidth=5)


    e_axis.set_ylabel('Eccentricity')
    w_axis.set_ylabel('Frequency [rad/day]')
    w_axis.set_xlabel('Age [Myr]')

#    if cmdline_config is not None:
#        w_axis.set_xlim(left=cmdline_config.disk_lifetime * 1e3)
#        e_axis.set_xlim(left=cmdline_config.disk_lifetime * 1e3)

    pyplot.figlegend()

def print_cmdline(cmdline_config):
    """Output to stdout the command line configuration for an evolution."""

    for key, value in vars(cmdline_config).items():
        print('\t* ' + repr(key) + ': ' + repr(value))

if __name__ == '__main__':
    cmdline_args = parse_command_line()
    to_file = dict(cmdline_args.to_file)
    with open(cmdline_args.pickle, 'rb') as pickle_file:
        evolution_index = 0
        while True:
            try:
                pickled_cmdline = pickle.load(pickle_file)
                pickled_info = pickle.load(pickle_file)
                pickled_evolution = pickle.load(pickle_file)
                print(80 * '*')
                print('* %-76.76s *' % ('Evolution #' + str(evolution_index)))
                print(80 * '*')
                print_cmdline(pickled_cmdline)
                print(pickled_info.format())
                plot_evolution(pickled_evolution,
                               pickled_cmdline,
                               cmdline_args.exclude_core_spin)
                if cmdline_args.to_file:
                    if str(evolution_index) in to_file:
                        print(pickled_evolution.format())
                        pyplot.savefig(to_file[str(evolution_index)])
                else:
                    pyplot.show()

                pyplot.clf()
                evolution_index += 1
            except EOFError:
                break
