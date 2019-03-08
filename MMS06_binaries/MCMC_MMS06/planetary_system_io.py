#!/usr/bin/python3 -u

#Legacy module name, too painful to change.
#pylint: disable=invalid-name
#pylint: enable=invalid-name

"""IO operations for planetary system parameters."""


import csv
from glob import glob
import re

from ast import literal_eval
import numpy
from astropy.units import Unit

from manual_exoplanet_data import data as manual_data

hatsouth_hostname = {'HATS610-018': 'HATS-18',
                     'HATS579-017': 'HATS-9'}

def add_io_arguments(parser, add_stellar_spin=False, add_mcmc=False):
    """Add arguments to a command line parser for finding the I/O files."""

    parser.add_argument(
        '--info-files',
        default=glob('inputs/versioned/*.info'),
        nargs='*',
        help='The names of the HATSouth info files containing planets that'
        ' should be used in the constraint. Default: %(default)s.'
    )
    parser.add_argument(
        '--nasa-planets',
        default='inputs/versioned/NASA_planets.csv',
        help='The name of the file containing information downloaded from '
        'the NASA exoplanet archive for planets that should be included in '
        'the constraint. Default: %(default)s.'
    )
    if add_stellar_spin:
        parser.add_argument(
            '--stellar-spin-file',
            default='inputs/versioned/stellar_spins.txt',
            help='The file from which to read the stellar spins. '
            'Default: %(default)s.'
        )
    if add_mcmc:
        parser.add_argument(
            '--mcmc-fname',
            default='outputs/unversioned/%(system)s.mcmc',
            help='A %%-substitution template giving the filename '
            'containincg the MCMC chain of system parameters produced by '
            'format_mcmc_samples.py. Default %(default)s.'
        )

#Type for objects holding data.
#pylint: disable=too-few-public-methods
class Structure:
    """An empty class used only to hold user defined attributes."""

    def __init__(self, **initial_attributes):
        """Create a class with (optionally) initial attributes."""

        for attribute_name, attribute_value in initial_attributes.items():
            setattr(self, attribute_name, attribute_value)

    def format(self, prefix=''):
        """Generate a tree-like representation of self."""

        result = ''
        for attr_name in dir(self):
            if attr_name[0] != '_':
                attribute = getattr(self, attr_name)
                if isinstance(attribute, Structure):
                    result += (prefix
                               +
                               '|-'
                               +
                               attr_name
                               +
                               '\n'
                               +
                               attribute.format(prefix + '| '))
                else:
                    result += (prefix
                               +
                               '|-'
                               +
                               attr_name
                               +
                               ': '
                               +
                               str(attribute)
                               +
                               '\n')
        return result
#pylint: enable=too-few-public-methods

class ArrayWithAttributes(numpy.ndarray):
    """Allow adding attributes to numpy arrays."""

def read_exoplanets_org(csv_filename):
    """
    Read a CSV file downloaded from exoplanets.org into a dictionary.

    Args:
        - csv_filename: The name of the comma separated file downloaded from
                        exoplanets.org.

    Returns: A dictionary with keys the column names in the input file and
             values being structures with attributes unit, giving the units
             specified in the input file and values, all the entries in the
             column converted to floating point values if appropriate.
    """

    string_data = csv.DictReader(open(csv_filename, 'r'))
    column_units = next(string_data)

    data = {k: [] for k in column_units.keys()}
    for data_row in string_data:
        for column, value in data_row.items():
            data[column].append(value)

    result = Structure()
    for k in data:
        try:
            setattr(result,
                    k,
                    numpy.array([
                        numpy.nan if v == '' else float(v) for v in data[k]
                    ]).view(ArrayWithAttributes))
        except ValueError:
            setattr(result,
                    k,
                    numpy.array(data[k]).view(ArrayWithAttributes))

    for k in column_units:
        getattr(result, k).unit = column_units[k]

    return result

def read_ages(nasa_planets,
              age_file_standard='inputs/versioned/getages.txt',
              age_file_manual_density='inputs/versioned/getages_nodensity.txt',
              manual_densities='inputs/versioned/age_variables_nodensity.txt'):
    """
    Complete the NASA exoplanet archive planets with age information.

    Args:
        - nasa_planets: The planets read from a CSV file downloaded from the
                        NASA exoplanet archive. On output, this gets updated
                        with the information from the various input files.
        - age_file_standard: The name of the file containing the derived
                             ages.
        - age_file_manual_density: The name of the file with ages derived
                                   from manually extracted densities.
        - manual_densities: The name of the file containing the manually
                            extracted densities themselves.

    Returns: None
    """

    def read_file(filename, columns):
        """
        Read one of the input files and update nasa_planets.

        Args:
            - filename: The name of the file to read.
            - columns: a dictionary of the quantities to read from the file
                       (keys) and the columns that contain them. The quantity
                       pl_hostname must be among the columns.

        Returns: None
        """

        hostname_list = list(nasa_planets.pl_hostname)
        num_systems = len(hostname_list)
        for quantity, column in columns.items():
            if not hasattr(nasa_planets, quantity):
                setattr(nasa_planets,
                        quantity,
                        numpy.full((num_systems,), numpy.nan))
        with open(filename, 'r') as input_file:
            for line in input_file:
                entries = line.split()
                host = entries[columns['pl_hostname']]
                system_index = 0
                while (
                        system_index < len(hostname_list)
                        and
                        (
                            not hostname_list[system_index].startswith(host)
                            or
                            (
                                len(hostname_list[system_index]) > len(host)
                                and
                                hostname_list[system_index][len(host)] != ' '
                            )
                        )
                ):
                    system_index += 1
                if system_index == len(hostname_list):
                    continue
                for quantity, column in columns.items():
                    if quantity != 'pl_hostname':
                        try:
                            entry_val = int(entries[column])
                        except ValueError:
                            entry_val = float(entries[column])
                        getattr(nasa_planets, quantity)[system_index] = (
                            entry_val
                        )

    age_file_columns = dict(pl_hostname=0,
                            st_mass=2,
                            st_masserr1=3,
                            st_rad=6,
                            st_raderr1=7,
                            st_age=10,
                            st_ageerr1=11,
                            st_lum=14,
                            st_lumerr1=15)

    density_file_columns = dict(pl_hostname=0,
                                st_dens=10,
                                st_denserr1=11,
                                st_denserr2=12)

    read_file(age_file_standard, age_file_columns)
    read_file(age_file_manual_density, age_file_columns)
    read_file(manual_densities, density_file_columns)

#Sufficient internal structure ensures simple enough code units.
#pylint: disable=too-many-locals
def read_nasa_planets(csv_filename,
                      eliminate=('SWEEPS-11',
                                 'HD 41004 B',
                                 'PSR J1719-1438',
                                 'K2-22'),
                      fill_missing=manual_data,
                      need_ages=True):
    """
    Read a CSV file downloaded from the NASA Exoplanet Archive to a dict.

    Args:
        csv_filename:    The name of the comma separated file downloaded from
            http://exoplanetarchive.ipac.caltech.edu.

    Returns:
        A dictionary with keys the column names in the input file and
        values being structures with attributes unit, giving the units
        specified in the input file and values, all the entries in the
        column converted to floating point values if appropriate.
    """

    def do_eliminate():
        """Eliminate the systems listed in eliminate."""

        system_names = [name.decode() for name in data[:, 1]]
        delete_indices = []
        for system in eliminate:
            if system in system_names:
                delete_indices.append(system_names.index(system))
        return numpy.delete(data, delete_indices, 0)

    def do_fill_missing(result):
        """Add the data from fill_missing to result."""

        system_names = list(result.pl_hostname)
        for fill_system in fill_missing:
            try:
                fill_index = system_names.index(fill_system['pl_hostname'])
            except ValueError:
                continue
            for quantity, value in fill_system.items():
                if hasattr(result, quantity):
                    getattr(result, quantity)[fill_index] = value

    data = numpy.genfromtxt(csv_filename,
                            delimiter=',',
                            comments='#',
                            dtype=None)

    data_columns = [col.decode() for col in data[0]]
    if eliminate:
        data = do_eliminate()

    result = Structure()
    column_name_list = []
    with open(csv_filename, 'r') as csv_file:
        for line in csv_file:
            if line[0] != '#':
                continue
            entries = line.strip().split()
            if len(entries) < 4 or entries[1] != 'COLUMN':
                continue
            column_name = entries[2].strip(':')
            column_index = data_columns.index(column_name)
            column_values = data[:, column_index][1:]
            if (
                    column_name in ['pl_hostname',
                                    'pl_name',
                                    'pl_discmethod',
                                    'pl_bmassprov',
                                    'st_optband',
                                    'rowupdate',
                                    'pl_letter',
                                    'pl_tsystemref',
                                    'pl_locale',
                                    'pl_facility',
                                    'pl_telescope',
                                    'pl_instrument',
                                    'pl_publ_date',
                                    'hd_name',
                                    'hip_name',
                                    'st_spstr',
                                    'st_metratio',
                                    'swasp_id']
                    or
                    column_name.endswith('_str')
                    or
                    column_name.endswith('link')
            ):
                column_values = [v.decode() for v in column_values]
            else:
                print('column_name: ' + repr(column_name))
                column_values = [numpy.nan if v == b'' else float(v)
                                 for v in data[:, column_index][1:]]
            setattr(
                result,
                column_name,
                numpy.array(column_values).view(ArrayWithAttributes)
            )
            column_name_list.append(column_name)

    if fill_missing:
        do_fill_missing(result)

    if need_ages:
        read_ages(result)

    for column_name in column_name_list:
        if column_name.endswith('err2'):
            column = getattr(result, column_name)
            nan_indices = numpy.isnan(column)
            column[nan_indices] = -getattr(result,
                                           column_name[:-1]+'1')[nan_indices]

    return result
#pylint: enable=too-many-locals

class FloatWithAttributes(float):
    """Allow adding attributes to floats."""

class StrWithAttributes(str):
    """Allow adding attributes to strings."""

class FloatWithErrors(float):
    """Floating point value with errors."""

    def __new__(cls, value, plus_error, minus_error=None):
        """Make a new float with plus_error and minus_error attributes."""

        result = super(FloatWithErrors, cls).__new__(cls, value)
        result.plus_error = plus_error
        result.minus_error = (plus_error
                              if minus_error is None else
                              minus_error)
        return result

    def __getnewargs__(self):
        #False positive, pylind does not see members created by __new__
        #pylint: disable=no-member
        return float(self), self.plus_error, self.minus_error
        #pylint: enable=no-member

    def __str__(self):
        """Show the value and errors."""

        return '%s +%s -%s' %(
            repr(float(self)),
            #False positive, pylind does not see members created by __new__
            #pylint: disable=no-member
            repr(self.plus_error),
            repr(self.minus_error)
            #pylint: enable=no-member
        )

def read_hatsouth_info(info_filename):
    """
    Read one of the info files with HATSouth system parameters.

    Args:
        info_filename:    The name of the info file to parse.

    Returns:
        A structure with all parameters contained in the info file,
        each of which is a floating point value but with attributes
        'plus_error' and 'minus_error' giving the error estimate in the
        positive and negative directions respectively and 'description'
        containing the description found in the info file.
    """

    result = Structure()
    with open(info_filename, 'r') as info:
        for line in info:
            if len(line) > 1 and line[0] != '%':
                definition = line.split('#')[0]
                name, value_str = definition.split('=')

                value_str = value_str.strip()
                if value_str[0] == '[' and value_str[-1] == ']':
                    value = literal_eval(value_str)
                else:
                    value_str = value_str.split()
                    assert len(value_str) >= 1

                    try:
                        value = float(value_str[0])
                    except ValueError:
                        value = value_str[0]

                #False positive, pylint does not see attributes
                #defined in __new__
                #pylint: disable=attribute-defined-outside-init
                if isinstance(value, float):
                    if len(value_str) != 1:
                        if len(value_str) == 4:
                            value *= Unit(value_str[3])
                        else:
                            value *= Unit('')

                        if value_str[1] == '+/-':
                            value.plus_error = float(value_str[2])
                            value.minus_error = value.plus_error
                        else:
                            value.plus_error = float(value_str[1])
                            value.minus_error = abs(float(value_str[2]))
                #pylint: enable=attribute-defined-outside-init


                setattr(result, name.strip(), value)

    return result

def read_stellar_spins(stellar_spin_fname):
    """Read the stellar spin period information from the given file."""

    spin_data = numpy.genfromtxt(stellar_spin_fname,
                                 names=True,
                                 dtype=None)
    return {
        system['Host'].decode(): FloatWithErrors(
            system['SpinPeriod'],
            system['SpinPeriodPlusError'],
            system['SpinPeriodMinusError']
        ) for system in spin_data
    }

def parse_hartman2016_stellar_params(filename):
    """Return table 7 from Hartman et. al. 2016 as a Structure with lists."""


    def read_by_word(skip_lines):
        """Iterate over the words in the file, since lines make no sense."""

        drop_mathrm = re.compile(r'\\(mathrm|rm)\{(?P<mathrm>.*?)\}')
        merge_colname = re.compile(r'  *_')
        skip_latex = str.maketrans('', '', '{}$\\')
        with open(filename, 'r') as param_file:
            for line in param_file:
                if line.startswith('References.'):
                    return
                if skip_lines > 0:
                    skip_lines -= 1
                else:
                    for word in merge_colname.sub(
                            '_',
                            drop_mathrm.sub(r'\g<mathrm>',
                                            line).translate(skip_latex)
                    ).split():
                        yield word

    def read_value_with_error(next_word, file_word):
        """
        Parse the various value formats to FloatWithErrors.

        Understands simple floating point values (e.g. 3.14159), floats with
        symmetric errors (e.g. '3.1415 +or- 0.0003') and floats with
        asymmetric errors (e.g. '${1.005}_{-0.027}^{+0.032}$')

        Args:
            - next_word:
                The next word in the file.

            - file_word:
                A file iterator which can be querried to get more words.

        Returns:
            - value:
                A FloatWithErrors instance contaning the parsed value.

            - next_word:
                The next word in the file.
        """

        asymmetric_error_rex = re.compile(
            r'(?P<value>[-+0123456789.]*)'
            r'_-(?P<minus_err>[0123456789.]*)'
            r'\^\+(?P<plus_err>[0123456789.]*)'
        )
        asymmetric_match = asymmetric_error_rex.match(next_word)
        if asymmetric_match is not None:
            return (
                FloatWithErrors(
                    float(asymmetric_match.group('value')),
                    float(asymmetric_match.group('plus_err')),
                    float(asymmetric_match.group('minus_err'))
                ),
                next(file_word)
            )
        if next_word == 'cdots':
            return (FloatWithErrors(float('nan'), float('nan')),
                    next(file_word))

        value = float(next_word)
        next_word = next(file_word)
        if next_word == '+or-':
            next_word = next(file_word)
            return (FloatWithErrors(value, float(next_word)),
                    next(file_word))

        return FloatWithErrors(value, float('nan')), next_word

    def read_reference_list(next_word, file_word):
        """Read a list of integers indicating refercences."""

        result = []

        if next_word == 'cdots':
            return result

        while next_word[-1] == ',':
            result.extend([int(r) for r in next_word[:-1].split(',')])
            next_word = next(file_word)
        result.extend([int(r) for r in next_word[:-1].split(',')])
        return result

    skip_lines = 2
    first_units = '(days)'
    surveys_3char = ['HAT', 'KEL', 'TrE', 'WAS']
    file_word = read_by_word(skip_lines)
    columns = []
    next_word = ''
    while next_word != first_units:
        next_word = next(file_word)
        if next_word != first_units:
            columns.append(next_word)

    result = Structure(**{c: [] for c in columns})

    while next_word[:3] not in surveys_3char:
        next_word = next(file_word)

    while True:
        assert next_word[:3] in surveys_3char
        try:
            for column_name in columns:
                if column_name == 'Planet':
                    #False positive, Planet member created out of __init__
                    #pylint: disable=no-member
                    result.Planet.append(next_word)
                    #pylint: enable=no-member
                    next_word = next(file_word)
                elif column_name == 'References':
                    #False positive, References member created out of __init__
                    #pylint: disable=no-member
                    result.References.append(
                        read_reference_list(next_word, file_word)
                    )
                    #pylint: enable=no-member
                    next_word = next(file_word)
                else:
                    value, next_word = read_value_with_error(next_word,
                                                             file_word)
                    getattr(result, column_name).append(value)
        except StopIteration:
            break

    return result

if __name__ == '__main__':
    stellar_params = parse_hartman2016_stellar_params(
        '../Scripts/inputs/versioned/'
        'stellar_parameters_Hartman_et_al_2016.txt'
    )

    print(stellar_params.format())
    print(
        'WASP-19 index: '
        +
        repr(
            next(
                #False positive, Planet member created out of __init__
                #pylint: disable=no-member
                (i for i in range(len(stellar_params.Planet))
                 if stellar_params.Planet[i].startswith('WASP-19'))
                #pylint: enable=no-member
            )
        )
    )
