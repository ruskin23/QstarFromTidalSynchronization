#!/usr/bin/env python3

import sys
import os
from pathlib import Path
from directories import directories

home_dir=str(Path.home())
path=directories(home_dir)
sys.path.append(path.poet_path+'/PythonPackage')
sys.path.append(path.poet_path+'/scripts')

from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as \
    orbital_evolution_library
from orbital_evolution.transformations import phase_lag

import pickle
import numpy
import scipy
import argparse
import matplotlib
from matplotlib import pyplot

from binary_evolution import Evolution

import time

import logging


def plot_evolution(evolved_binary):

    wsun = 0.24795522138
    evolution=evolved_binary.get_evolution()

    wenv_secondary = (evolution.secondary_envelope_angmom / evolved_binary.secondary.envelope_inertia(evolution.age)) / wsun
    wcore_secondary = (evolution.secondary_core_angmom / evolved_binary.secondary.core_inertia(evolution.age)) / wsun
    wenv_primary = (evolution.primary_envelope_angmom / evolved_binary.primary.envelope_inertia(evolution.age)) / wsun
    wcore_primary = (evolution.primary_core_angmom / evolved_binary.primary.core_inertia(evolution.age)) / wsun
    orbitalfrequncy = evolved_binary.orbital_frequency(evolution.semimajor) / wsun

    pyplot.semilogx(evolution.age, wenv_primary, color="b", label='Primary Star Envelope')
    pyplot.semilogx(evolution.age, wenv_secondary, color="r", label='Secondary Star Envelope')
    pyplot.semilogx(evolution.age, wcore_primary, color="b", linestyle='--', label='Primary Star Core')
    pyplot.semilogx(evolution.age, wcore_secondary, color="r",linestyle='--', label='Secondary Star Core')

    pyplot.semilogx(evolution.age, orbitalfrequncy, "-k", label='Orbital Frequency')
    pyplot.legend(loc='upper right')
    pyplot.ylabel('Spin Freuqncy')
    pyplot.xlabel('age')
    pyplot.ylim(top=100)
    pyplot.ylim(bottom=-20)
    pyplot.show()

def get_initial_conditions(evolution):

    results = evolution.calculate_intial_conditions()
    return results

def get_evolution(evolution):

    evolved_binary=evolution.evolve_binary()
    final_state=evolved_binary.final_state()
    spin = (2.0*numpy.pi*evolved_binary.primary.envelope_inertia(final_state.age)/
                final_state.primary_envelope_angmom)
    final_orbital_period=evolved_binary.orbital_period(final_state.semimajor)
    final_eccentricity=final_state.eccentricity

    print('Final Orbital Period = {} \nFinal Eccentricity = {} \n Final Spin ={}'.format(final_orbital_period,final_eccentricity,spin))
    return spin,final_orbital_period,final_eccentricity
    # plot_evolution(evolved_binary)



if __name__=='__main__':

    #clear log file. this is best way I found
    with open('spin_calculation.log', 'w'): pass


    parser = argparse.ArgumentParser()
    parser.add_argument('-s',
                        dest='system',
                        default=None,
                        help='pick system from a file')

    args = parser.parse_args()

    logging.basicConfig(filename='spin_calculation.log',level=logging.INFO)

    serialized_dir = path.poet_path +  "/stellar_evolution_interpolators"
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    eccentricity_path=os.path.join(path.poet_path,'eccentricity_expansion_coef.txt').encode('ascii')

    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        eccentricity_path
    )

    parameters=dict()


    if args.system is not None:
        parameters['system']=args.system
        with open('SpinlogQCatalog_el0.4.txt','r') as f:
            for lines in f:
                x=lines.split()
                if x[0]==args.system:
                    primary_mass=float(x[15])
                    massratio=float(x[14])
                    secondary_mass=massratio*primary_mass
                    age=float(x[16])
                    feh=float(x[17])
                    spin_period=float(x[12])
                    orbital_period=float(x[6])
                    eccentricity=float(x[8])
                    Wdisk=4.1
                    logQ=12.0

    parameters['primary_mass']=primary_mass
    parameters['secondary_mass']=secondary_mass
    parameters['age']=age
    parameters['feh']=feh
    parameters['orbital_period']= orbital_period
    parameters['eccentricity']=eccentricity
    parameters['spin_period']=spin_period
    parameters['logQ']=logQ
    parameters['Wdisk']=Wdisk

    parameters['dissipation']=True
    parameters['disk_dissipation_age']=5e-3
    parameters['wind']=True
    parameters['wind_saturation_frequency']=2.54
    parameters['diff_rot_coupling_timescale']=5e-3
    parameters['wind_strength']=0.17
    parameters['print_cfile']=False
    parameters['evolution_max_time_step']=1e-3
    parameters['evolution_precision']=1e-6
    parameters['inclination']=0.0
    parameters['spin_frequency_breaks']=None
    parameters['spin_frequency_powers']=numpy.array([0.0])
    parameters['tidal_frequency_breaks']=None
    parameters['tidal_frequency_powers']=numpy.array([0.0])

    print(parameters)

    evolution=Evolution(interpolator,parameters)
    results =  evolution.calculate_intial_conditions()
    spin=results['spin']
    print('Spin Period = ',spin)









