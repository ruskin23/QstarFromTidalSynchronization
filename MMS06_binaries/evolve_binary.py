#!/usr/bin/env python3

from argparse import ArgumentParser

def parse_command_line():
    """Parse the command line arguments into object with attributes."""

    parser = ArgumentParser(
        description='Calculate orbital evolution of one of the MMS06 binaries.'
    )
    parser.add_argument(
        'binary',
        help='The info file for the binary to evolve.'
    )
    parser.add_argument(
        '--lgQ',
        type=float,
        default='6.0',
        help='The value of log10(Q*) to assume for both the primary and the '
        'secondary. Default: %(default)s.'
    )
    return parser.parse_args()

if __name__ == '__main__':
    cmdline_args = parse_command_line()
