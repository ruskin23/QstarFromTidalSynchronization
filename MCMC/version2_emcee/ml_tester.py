#!/usr/bin/env python3

from mcmc_sampler import sampler

import numpy

import os
import sys

from pathlib import Path
from directories import directories


home_dir=str(Path.home())
path=directories(home_dir)
sys.path.append(path.poet_path+'/PythonPackage')
sys.path.append(path.poet_path+'/scripts')
sys.path.append(path.binary_directory)

from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as \
    orbital_evolution_library

import multiprocessing
import logging


def initialize_interpolator():

    serialized_dir = path.poet_path +  "/stellar_evolution_interpolators"
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    eccentricity_path=os.path.join(path.poet_path,'eccentricity_expansion_coef_O400.sqlite').encode('ascii')

    orbital_evolution_library.prepare_eccentricity_expansion(
        eccentricity_path,
        1e-4,
        True,
        True
    )

    return interpolator

def sample_parameters(logger, kic):

    unit_cube = numpy.random.rand(9)
    print(unit_cube)
    logger.info('Begin Conversion for unit cube = {}'.format(unit_cube))
    transform = sampler(kic, unit_cube)
    parameters = transform()
    logger.info('Converted Parameters = {}'.format(parameters))

def test_sample(kic):

    unit_cube = numpy.random.rand(9)
    print(unit_cube)
    transform = sampler(kic, unit_cube)
    parameters = transform()

def worker(num, kic):

    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
        handler.close()

    logger = logging.getLogger('Process {}'.format(num))
    logger.setLevel(logging.INFO)
    file_hander = logging.FileHandler('process_{}.log'.format(multiprocessing.current_process().pid))
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_hander.setFormatter(formatter)
    logger.addHandler(file_hander)

    logger.info('Starting process {}'.format(num))
    sample_parameters(logger, kic)
    logger.info('Ending Process {}'.format(num))

if __name__ == '__main__':

    # interpolator = initialize_interpolator()
    kic = '10960995'

    for i in range(5):
        test_sample(kic)

    # p1 = multiprocessing.Process(target=worker, args=(1, kic,))
    # p1.start()
    # p1.join()

    # p2 = multiprocessing.Process(target=worker, args=(2, kic,))
    # p2.start()
    # p2.join()
