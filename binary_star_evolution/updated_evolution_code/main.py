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
                    # logQ=7.0

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
    parameters['tidal_frequency_breaks']=None
    parameters['tidal_frequency_powers']=numpy.array([0.0])

    Q_values=[5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,11.0,12]
    
    with open(f'spin_vs_logQ_{breaks}.txt','w') as fr:
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
    # get_evolution(evolution)




    # massratio=0.663
    # with open('results.txt','w') as ff:
    #     ff.write('porb_initial\te_initial\tWdisk\tlogQ\tprimary_mass\tage\tfeh\tporb_final\te_final\tspin\tfail\n')
    #     with open('temp.out','r') as f:
    #         for lines in f:
    #             p=lines.split()
    #             print(float(p[9][:-1]))
    #             parameters['system']='1'
    #             parameters['primary_mass']=float(p[9][:-1])
    #             parameters['secondary_mass']=float(p[9][:-1])*massratio
    #             parameters['age']=float(p[11][:-1])
    #             parameters['feh']=float(p[13][:-1])
    #             parameters['orbital_period']= float(p[1][:-1])
    #             parameters['eccentricity']=float(p[3][:-1])
    #             parameters['logQ']=float(p[7][:-1])
    #             parameters['Wdisk']=float(p[5][:-1])


    #             for item,value in parameters.items():
    #                 print('{} = {}'.format(item,value))
            
    #             evolution=Evolution(interpolator,parameters)
    #             evolved_binary=evolution.evolve_binary()
    #             final_state=evolved_binary.final_state()
    #             if final_state.age!=parameters['age']:
    #                 spin,final_orbital_period,final_eccentricity=scipy.nan,scipy.nan,scipy.nan
    #             else:
    #                 spin = (2.0*numpy.pi*evolved_binary.primary.envelope_inertia(final_state.age)/
    #                         final_state.primary_envelope_angmom)
    #                 final_orbital_period=evolved_binary.orbital_period(final_state.semimajor)
    #                 final_eccentricity=final_state.eccentricity

    #             ff.write(repr(float(p[1][:-1]))+'\t'+
    #                      repr(float(p[3][:-1]))+'\t'+
    #                      repr(float(p[5][:-1]))+'\t'+
    #                      repr(float(p[7][:-1]))+'\t'+
    #                      repr(float(p[9][:-1]))+'\t'+
    #                      repr(float(p[11][:-1]))+'\t'+
    #                      repr(float(p[13][:-1]))+'\t'
    #                      )
    #             ff.write(repr(final_orbital_period)+'\t'+repr(final_eccentricity)+'\t'+repr(spin)+'\t'+repr(numpy.isnan(spin))+'\n')








    # with open('parameters.txt','r') as f:
    #     next(f)
    #     for lines in f:
    #         x=lines.split()
    #         orbital_period=float(x[0])
    #         feh=float(x[1])
    #         eccentricity=float(x[2])
    #         Wdisk=float(x[3])
    #         logQ=float(x[4])
    #         primary_mass=float(x[5])
    #         age=float(x[6])
    #         Wdisk=float(x[7])
    #         secondary_mass=primary_mass*0.606
    #         spin_period=14.581
                
    # with open('/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/version1/SAVED_CHAINS/ganymede/MCMC_1/accepted_parameters_1.txt','r') as f:
    #     next(f)
    #     for lines in f:
    #         x=lines.split()
    #         if x[0]=='1':
    #             orbital_period=float(x[11])
    #             eccentricity=float(x[12])
    #             Wdisk=float(x[3])
    #             logQ=float(x[4])
    #             primary_mass=float(x[5])
    #             age=float(x[6])
    #             feh=float(x[7])
    #             secondary_mass=float(x[8])
    #             spin_period=10.778
    #             break








    """
    ages=[0.1,0.5,1.0]

    for age in ages:
        parameters['age']=age

        for key,value in parameters.items():
            print("{} = {}".format(key,value))

        for logQ1 in q:
            result=dict()
            if args.breaks:

                alpha=-1.0

                logQMax=4.0
                phase_lagMax=phase_lag(logQMax)
                phase_lag1=phase_lag(logQ1)
                omegaref1=2*numpy.pi

                omegaref=omegaref1*((phase_lagMax/phase_lag1)**(1/alpha))

                TidalFrequencyBreaks=numpy.array([omegaref])
                TidalFrequencyPowers=numpy.array([0,alpha])
                parameters['breaks']=True
                parameters['logQ']=logQMax
            else:
                TidalFrequencyBreaks=None
                TidalFrequencyPowers=numpy.array([0.0])
                parameters['breaks']=False
                parameters['logQ']=logQ1

            parameters['tidal_frequency_breaks']=TidalFrequencyBreaks
            parameters['tidal_frequency_powers']=TidalFrequencyPowers

            print(TidalFrequencyBreaks)

            evolve=Evolution(interpolator,parameters)
            final_values,frequencies=evolve()

        observed_spin_frquency=(2*numpy.pi/PspinCurrent)/wsun
        calculated_spin_frequencies=frequencies['wenv_primary']
        spin_diff=calculated_spin_frequencies-observed_spin_frquency

        print('Observed Spin = ',observed_spin_frquency)
        print('Caclulated Spin Frequency = ', calculated_spin_frequencies)
        print('Spin Diff = ', spin_diff)


        zero_crossing=numpy.where(numpy.diff(numpy.sign(spin_diff)))[0]

        print('Zero Crossing = ', zero_crossing)

        print(calculated_spin_frequencies[zero_crossing])
        print(frequencies['age'][zero_crossing])


    #checking phase lag bechaviour with tidal frequency:
    #phaseLag=[]
    #for omega in result['tidal_frequency']:
    #    phaseLag.append(phase_lagMax*(omegaref/omega))

    #pyplot.plot(result['tidal_frequency'],phaseLag)
    #pyplot.show()

    if parameters['plot']==True:


        pyplot.axhline(y=(2*numpy.pi/parameters['PspinCurrent'])/wsun,linestyle=':',label='PSpinCurrent',color='r')
        pyplot.axhline(y=2*numpy.pi/parameters['PorbCurrent']/wsun,linestyle=':',label='PorbCurrent',color='g')
        pyplot.legend()
        pyplot.show()
    """

