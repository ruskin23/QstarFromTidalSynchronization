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

    parser = argparse.ArgumentParser()
    parser.add_argument('-s',
                        dest='system',
                        default=None,
                        help='pick system from a file')

    parser.add_argument('-b',
                        dest='breaks',
                        default=None,
                        help='decide if breaks needed or not')

    args = parser.parse_args()
    breaks=float(args.breaks)


    serialized_dir = path.poet_path +  "/stellar_evolution_interpolators"
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    eccentricity_path=os.path.join(path.poet_path,'eccentricity_expansion_coef.txt').encode('ascii')

    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        eccentricity_path
    )

    parameters=dict()


    parameters['system']=args.system
    with open(path.current_directory+'/SpinlogQCatalog_el0.4.txt','r') as f:
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

    parameters['primary_mass']=primary_mass
    parameters['secondary_mass']=secondary_mass
    parameters['age']=age
    parameters['feh']=feh
    parameters['orbital_period']= orbital_period
    parameters['eccentricity']=eccentricity
    parameters['spin_period']=spin_period
    # parameters['logQ']=logQ
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
    parameters['breaks']=True
    parameters['spin_frequency_breaks']=None
    parameters['spin_frequency_powers']=numpy.array([0.0])

    Q_values=[5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,11.0,12]
    
    with open(f'{path.results_directory}/System_{args.system}/spin_vs_logQ_{breaks}.txt','w') as fr:
        fr.write('logQ\tspin_constant_Q\tspin_freq_Q\tp_initial\te_initial\tp_final\te_final\tdelta_p\tdelta_e\n')
        for q in Q_values:
            print(f'Solving for logQ = {q}')
            alpha=breaks
            logQMax=4.0
            logQ1=q
            phase_lagMax=phase_lag(logQMax)
            phase_lag1=phase_lag(logQ1)
            omegaref1=2*numpy.pi
            omegaref=omegaref1*((phase_lagMax/phase_lag1)**(1/alpha))
            print(omegaref)
            print(f'observed spin period={spin_period}')
            parameters['logQ']=logQ1
            parameters['tidal_frequency_breaks']=numpy.array([omegaref])
            parameters['tidal_frequency_powers']=numpy.array([0,alpha])
            evolution=Evolution(interpolator,parameters)
            results =  evolution.calculate_intial_conditions()
            fr.write(repr(q)+'\t'+
                    repr(results['spin'])+'\t'+
                    repr(results['p_initial'])+'\t'+
                    repr(results['e_initial'])+'\t'+
                    repr(results['p_final'])+'\t'+
                    repr(results['e_final'])+'\t'+
                    repr(results['delta_p'])+'\t'+
                    repr(results['delta_e'])+'\n'
                    )
