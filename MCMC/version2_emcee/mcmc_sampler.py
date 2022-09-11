# from distutils.command.config import config
from distutils.command.config import config
from configargparse import ArgumentParser
from gettext import Catalog
import logging
import sys
import os
import random
import numpy
import scipy

from pathlib import Path
from directories import directories


home_dir=str(Path.home())
path=directories(home_dir)
sys.path.append(path.poet_path+'/PythonPackage')
sys.path.append(path.poet_path+'/scripts')
sys.path.append(path.binary_directory)


from stellar_evolution.library_interface import library
from create_objects import BinaryObjects
from prior_transform_class import prior_transform
from nested_solver import InitialConditionSolver
from initial_secondary_angmom import IntialSecondaryAngmom

from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as \
    orbital_evolution_library
from orbital_evolution.transformations import phase_lag
from multiprocessing import Pool
from Logger import setup_process


import emcee
from  hdf5_backend import HDFBackend

_logger=logging.getLogger(__name__)


def cmdline_parser():

    p = ArgumentParser(default_config_files=['config.txt'])
    p.add_argument('--system')
    p.add_argument('--num_parallel_processes',type=int,default=16,help='number of parallel processes')
    # p.add_argument('--fname_datetime_format',default='%Y%m%d%H%M%S')
    p.add_argument('--fname_datetime_format',default='%Y%m%d%H%M%S')
    p.add_argument('--logging_datetime_format',default=None)
    p.add_argument('--std_out_err_fname',default='%(system)s_%(now)s_%(pid)d.err')
    p.add_argument('--logging_message_format',default=None)
    p.add_argument('--logging_fname',default='%(system)s_%(now)s_%(pid)d.log')

    return p.parse_args()


class UnchunkedPool:
    """Disable chunking in Pool.map."""

    def __init__(self, pool):
        """Wrap around the given pool's map."""

        self._pool = pool

    def map(self, *args, **kwargs):
        """Delegate everything to parent, but set chunksize to 1."""

        return self._pool.map(*args, **kwargs, chunksize=1)


class sampler:

    def __init__(self,
                system_number,
                uniform_variables):

        self.system_number=system_number
        self.uniform_variable=uniform_variables

        self.params=dict()

    def param_conversion(self,
                        sampled_params):

        sum_mass=sampled_params[1]
        mass_ratio=sampled_params[2]

        primary_mass=(1/(1+mass_ratio))*sum_mass
        secondary_mass=primary_mass*mass_ratio

        return [primary_mass,secondary_mass]


    def sampled_from_coupled(self):

        params=prior_transform(self.system_number)
        sampled_params=params.get_samples(self.uniform_variable)
        masses=self.param_conversion(sampled_params)

        self.params['primary_mass'] = masses[0]
        self.params['secondary_mass'] = masses[1]
        z=sampled_params[0]
        self.params['feh'] = library.feh_from_z(z)
        self.params['age'] = (10**(sampled_params[3]))/1e9
        self.params['eccentricity'] = sampled_params[4]



    def fixed_params(self):

        self.params['dissipation']=True
        self.params['disk_dissipation_age']=5e-3
        self.params['wind']=True
        self.params['wind_saturation_frequency']=2.54
        self.params['diff_rot_coupling_timescale']=5e-3
        self.params['wind_strength']=0.17

    def uniformly_sampled(self):

        phase_lag_sampled = phase_lag( ( 7*self.uniform_variable[5] ) + 5 )
        alpha = ( 10*self.uniform_variable[6] ) - 5
        lnP = ( numpy.log(50)-numpy.log(0.5) )*self.uniform_variable[7] + numpy.log(0.5)

        omegaref = numpy.exp(numpy.log(2*numpy.pi) - lnP)
        omegamin =(2*numpy.pi)/50

        if alpha<0:
            tidal_power = numpy.array([0,alpha])
            self.params['tidal_frequency_breaks'] = numpy.array([omegaref])
            phase_lag_max=phase_lag_sampled
        else:
            tidal_power = numpy.array([0,alpha,0])
            self.params['tidal_frequency_breaks'] = numpy.array([omegamin,omegaref])
            phase_lag_max = phase_lag_sampled*((omegamin/omegaref)**alpha)

        self.params['phase_lag_max']=phase_lag_max
        self.params['tidal_frequency_powers']=tidal_power
        self.params['spin_frequency_breaks']=None
        self.params['spin_frequency_powers']=numpy.array([0.0])

        disk_N = (1.0/1.4) - (1.0/14)
        self.params['Wdisk'] = 2 * numpy.pi * ( ( disk_N  * self.uniform_variable[8] )  +  1.0/14.0 )
        _logger.info('Wdisk = {!r}'.format(self.params['Wdisk']))

        return alpha,omegaref

    def get_orbital_period(self):
        with open(path.current_directory+'/catalog/filtering/nominal_value_catalog_Iconv_cutoff.txt','r') as f:
            for lines in f:
                x=lines.split()
                if self.system_number==x[1]:
                    self.params['orbital_period']=float(x[2])
                    break


    def __call__(self):

        alpha,omegaref=self.uniformly_sampled()
        self.fixed_params()
        self.sampled_from_coupled()
        self.get_orbital_period()
        self.params['function']='nested'
        self.params['method']='brentq'

        return self.params,alpha,omegaref


def log_probablity(unit_cube_values,interpolator,system_number,observed_spin):

    invalid_values = tuple( -numpy.inf if i == 0 else numpy.nan
                        for i in range(len(unit_cube_values)+1)
                        )

    if unit_cube_values.min() < 0 or unit_cube_values.max() > 1:
        _logger.warning(
            'At least one proposed unit cube value is outside the range(0, 1): '
            '%s',
            repr(unit_cube_values)
        )
        return invalid_values


    _logger.info('Begin Conversion using %s',repr(unit_cube_values))

    sampled_params=sampler(system_number,unit_cube_values)
    parameter_set,alpha,omegaref=sampled_params()

    _quantities=['primary_mass',
                'secondary_mass',
                'feh',
                'age',
                'Wdisk',
                'phase_lag_max',
                'tidal_frequency_breaks',
                'tidal_frequency_powers']

    _logger.info('\nFor debugging')
    for key in _quantities:
        _logger.info('parameters["{}"]={}'.format(key,parameter_set[key]))


    if parameter_set['feh'] < -1.014 or parameter_set['feh'] > 0.53:
        _logger.warning('feh value = %s is out of POET range (-1.014,0.53)',repr(parameter_set['feh']))
        return invalid_values

    if parameter_set['primary_mass'] < 0.4 or parameter_set['primary_mass'] > 1.2:
        _logger.warning('primary mass value = %s is out of POET range (0.4,1.2)',repr(parameter_set['primary_mass']))
        return invalid_values

    if parameter_set['secondary_mass'] < 0.4 or parameter_set['secondary_mass'] > 1.2:
        _logger.warning('secondary mass value = %s is out of POET range (0.4,1.2)',repr(parameter_set['secondary_mass']))
        return invalid_values


    for _quantity in ['primary','secondary']:
        age_max=interpolator('radius',
                             parameter_set[_quantity+'_mass'],
                             parameter_set['feh']
                             ).max_age
        if parameter_set['age']>age_max:
            _logger.info('Sampled age is greater than max age allowed from interpolator. Setting sampled_age=age_max-1e-4')
            parameter_set['age']=age_max-1e-4


    _logger.info('\nSampled Parameters:')
    for key,values in parameter_set.items():
        _logger.info(f'{key}\t{values}\n')

    stars=BinaryObjects(interpolator,parameter_set)
    primary=stars.create_star(parameter_set['primary_mass'])
    secondary=stars.create_star(parameter_set['secondary_mass'])

    angmom=IntialSecondaryAngmom(interpolator,parameter_set)

    initial_conditions=InitialConditionSolver(interpolator,parameter_set,secondary_angmom=angmom())

    spin=initial_conditions(primary,secondary)

    log_likelihood=scipy.stats.norm(observed_spin['value'],observed_spin['sigma']).logpdf(spin)

    _logger.info('log_likelihood={!r}'.format(log_likelihood))

    p_names=['primary_mass',
             'secondary_mass',
             'feh',
             'age',
             'eccentricity',
             'Wdisk',
             'phase_lag_max']
    parameters=tuple(parameter_set[param_name] for param_name in p_names) + (alpha,omegaref)

    return ((log_likelihood,) + parameters)



if __name__ == '__main__':

    print('main program initiated')

    config=cmdline_parser()
    setup_process(config)

    print('Initializing interpolator')
    serialized_dir =  path.poet_path+"/stellar_evolution_interpolators"
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')


    eccentricity_path=os.path.join(path.scratch_directory,'eccentricity_expansion_coef_O400.sqlite').encode('ascii')

    orbital_evolution_library.prepare_eccentricity_expansion(
        eccentricity_path,
        1e-4,
        True,
        True
    )

    print('Interpolator initialized')


    system_number=config.system
    observed_spin=dict()
    with open(path.current_directory+'/catalog/filtering/Lurie_binaries_with_p1.txt','r') as f:
        for lines in f:
            x=lines.split()
            KIC=x[0]
            if KIC==system_number:
                observed_spin['value']=float(x[5])
                observed_spin['sigma']=abs(float(x[5])-float(x[6]))
                break

    #Initialising the state from 8-parameter set
    h5_filename = path.scratch_directory+'/sampling_output/h5_files/8Dfiles/system_'+system_number+'.h5'
    reader = emcee.backends.HDFBackend(h5_filename, read_only=True)
    last_state = reader.get_last_sample().coords
    wdisk_initial_state = numpy.random.rand(64,1)
    initial_state = numpy.array([numpy.append(last_state[i],wdisk_initial_state[i]) for i in range(64)])

    #h5 file to save and continue sampling
    backend_reader = HDFBackend(path.scratch_directory+'/sampling_output/h5_files'+'/system_'+system_number+'.h5')

    _logger.info('\nInitial State generated: ')
    _logger.info(initial_state)

    #Initiate blobs
    parameters=['primary_mass','secondary_mass', 'feh','age','eccentricity','Wdisk','phase_lag_max','alpha','break_period']
    blobs_dtype = [(name, float) for name in parameters]
    blobs_dtype = numpy.dtype(blobs_dtype)

    nwalkers = 64
    ndim = 9

    with Pool(
            config.num_parallel_processes,
            initializer=setup_process,
            initargs=[config],
            maxtasksperchild=1
    ) as workers:

        sampler_emcee=emcee.EnsembleSampler(nwalkers,
                                            ndim,
                                            log_probablity,
                                            args=(interpolator,system_number,observed_spin),
                                            blobs_dtype=blobs_dtype,
                                            backend=backend_reader,
                                            pool=UnchunkedPool(workers)
                                            )

        sampler_emcee.run_mcmc(None,nsteps=1000,progress=False)


