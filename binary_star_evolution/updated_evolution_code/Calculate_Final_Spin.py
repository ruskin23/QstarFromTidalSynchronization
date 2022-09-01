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


import numpy
import argparse
from matplotlib import pyplot

from binary_evolution import Evolution

import logging



def get_parameters_logfile(logfilename):

    _simple_quantities=['primary_mass',
                        'secondary_mass',
                        'feh',
                        'age',
                        'Wdisk',
                        'orbital_period',
                        'eccentricity',
                        'phase_lag_max']

    parameters=dict()

    with open('logfile/files/'+logfilename,'r') as f:
        for i,lines in enumerate(f):

            if i>6 and i<42 and i%2==1:

                x=lines.split()

                if x[0] in _simple_quantities:
                    parameters[x[0]]=float(x[1])

                if x[0]=='tidal_frequency_breaks':
                    x=lines.split()
                    if len(x)==2:
                        value=float(x[1][1:-1])
                    elif len(x)==3:
                        a=float(x[1][1:])
                        b=float(x[-1][:-1])
                        value=numpy.array([a,b])
                    parameters[x[0]]=numpy.array(value)

                if x[0]=='tidal_frequency_powers':
                    if len(x)==4:
                        a=float(x[2])
                        b=float(x[3][0:-1])
                        value=numpy.array([a,b])
                    elif len(x)==5:
                        a=float(x[1][1:])
                        b=float(x[2])
                        c=float(x[3])
                        value=numpy.array([a,b,c])
                    parameters[x[0]]=numpy.array(value)
        
    return parameters


def plot_evolution(evolved_binary):

    wsun = 0.24795522138
    evolution=evolved_binary.get_evolution()

    wenv_secondary = evolution.secondary_envelope_angmom / evolved_binary.secondary.envelope_inertia(evolution.age)
    wcore_secondary = evolution.secondary_core_angmom / evolved_binary.secondary.core_inertia(evolution.age)
    wenv_primary = evolution.primary_envelope_angmom / evolved_binary.primary.envelope_inertia(evolution.age)
    wcore_primary = evolution.primary_core_angmom / evolved_binary.primary.core_inertia(evolution.age)
    orbitalfrequncy = evolved_binary.orbital_frequency(evolution.semimajor)

    pyplot.semilogx(evolution.age, 2*numpy.pi/wenv_primary, color="b", label='Primary Star Envelope')
    pyplot.semilogx(evolution.age, 2*numpy.pi/wenv_secondary, color="r", label='Secondary Star Envelope')
    # pyplot.semilogx(evolution.age, wcore_primary, color="b", linestyle='--', label='Primary Star Core')
    # pyplot.semilogx(evolution.age, wcore_secondary, color="r",linestyle='--', label='Secondary Star Core')

    pyplot.semilogx(evolution.age, 2*numpy.pi/orbitalfrequncy, "-k", label='Orbital Frequency')
    pyplot.legend(loc='upper right')
    pyplot.ylabel('Spin Freuqncy')
    pyplot.xlabel('age')
    # pyplot.ylim(top=100)
    # pyplot.ylim(bottom=-20)
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

    parser = argparse.ArgumentParser()
    parser.add_argument('-s',
                        dest='system',
                        default=None,
                        help='pick system from a file')

    parser.add_argument('-f')
    parser.add_argument('-m')


    args = parser.parse_args()

    with open(path.scratch_directory+'/spin_calculation_{}_{}.log'.format(args.f,args.m),'w'):
        pass
    logging.basicConfig(filename=path.scratch_directory+'/spin_calculation_{}_{}.log'.format(args.f,args.m),
                        level=logging.INFO)

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

    parameters=dict()

    parameters['function']=args.f
    parameters['method']=args.m

    parameters['dissipation']=True
    parameters['wind']=True
    parameters['disk_dissipation_age']=5e-3
    parameters['wind_saturation_frequency']=2.54
    parameters['diff_rot_coupling_timescale']=5e-3
    parameters['wind_strength']=0.17
    parameters['evolution_max_time_step']=1e-3
    parameters['evolution_precision']=1e-6
    parameters['inclination']=0.0
    parameters['spin_frequency_breaks']=None
    parameters['spin_frequency_powers']=numpy.array([0.0])


    parameters["primary_mass"]=0.9241915970492054
    parameters["secondary_mass"]=0.8593663624245791
    parameters["feh"]=0.46458608151644215
    parameters["age"]=6.289791713929844
    parameters["Wdisk"]=1.1954168402393623
    parameters["phase_lag_max"]=6.251364607650211e-09
    parameters["tidal_frequency_breaks"]=numpy.array([0.12566371,7.02469755])
    parameters["tidal_frequency_powers"]=numpy.array([0.0,1.51411758, 0.0])

    parameters['orbital_period']= 11.653641112734459
    parameters['eccentricity']= 0.4406876762698404


    age_max_m1=interpolator('radius', parameters['primary_mass'], parameters['feh']).max_age
    print('max_age_allowed {} target_age {}'.format(age_max_m1,parameters['age']))

    print(parameters)

    evolution=Evolution(interpolator,parameters)
    # evolved_binary=evolution.evolve_binary()
    # plot_evolution(evolved_binary)
    spin =  evolution.calculate_intial_conditions()
    print('Spin Period = ',spin)