#!/usr/bin/env python3

import sys
import os
from pathlib import Path
from tkinter.messagebox import NO
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
    parser.add_argument('-l',
                        dest='logfile',
                        default=None,
                        help='select logfile to get parameters')

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

    eccentricity_path=os.path.join(path.scratch_directory,'eccentricity_expansion_coef_O400.sqlite').encode('ascii')

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
    parameters['disk_dissipation_age']=5e-3
    parameters['wind']=True
    parameters['wind_saturation_frequency']=2.54
    parameters['diff_rot_coupling_timescale']=5e-3
    parameters['wind_strength']=0.17
    parameters['print_cfile']=False
    parameters['evolution_max_time_step']=1e-2
    parameters['evolution_precision']=1e-5
    parameters['inclination']=0.0
    parameters['spin_frequency_breaks']=None
    parameters['spin_frequency_powers']=numpy.array([0.0])



    if args.logfile is not None:
        _simple_quantities=['primary_mass',
                    'secondary_mass',
                    'feh',
                    'age',
                    'Wdisk',
                    'orbital_period',
                    'eccentricity',
                    'phase_lag_max']


        with open('logfile/files/'+args.logfile,'r') as f:
            for i,lines in enumerate(f):

                if i>6 and i<42 and i%2==1:

                    x=lines.split()

                    if x[0] in _simple_quantities:
                        parameters[x[0]]=float(x[1])

                    if x[0]=='tidal_frequency_breaks':
                        x=lines.split()
                        if len(x)==2:
                            parameters[x[0]]=numpy.atleast_1d(float(x[1][1:-1]))
                        elif len(x)==3:
                            a=float(x[1][1:])
                            b=float(x[-1][:-1])
                            parameters[x[0]]=numpy.array([a,b])
                        elif len(x)==4:
                            a=float(x[1][1:])
                            b=float(x[-2])
                            parameters[x[0]]=numpy.array([a,b])


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


    age_max_m1=interpolator('radius', parameters['primary_mass'], parameters['feh']).max_age
    print(age_max_m1,parameters['age'])
    #Bad Solution for both dp and de
    # parameters['primary_mass']=1.0780738620002064
    # parameters['secondary_mass']=1.0256723503503944
    # parameters['age']=3.452331590976395
    # parameters['feh']= -0.18318213164525268
    # parameters['orbital_period']= 4.14387936268

    # parameters['eccentricity']=0.005117792869355689
    # parameters['spin_period']=10
    # parameters['Wdisk']= 1.4460263992995446

    # parameters['phase_lag_max']=1.2336514749258355e-08
    # parameters['spin_frequency_breaks']=None
    # parameters['spin_frequency_powers']=numpy.array([0.0])
    # parameters['tidal_frequency_breaks']=numpy.array([0.12566371,0.30692368])
    # parameters['tidal_frequency_powers']=numpy.array([ 0.0,0.19434417,0.0])



    #Final eccentricity<1e-8
    # parameters['primary_mass']=1.1993719870830584
    # parameters['secondary_mass']=1.1462126964961248
    # parameters['age']= 1.2135404663745835
    # parameters['feh']=  -0.03896229835997274
    # parameters['orbital_period']= 4.14387936268

    # parameters['eccentricity']=  0.0035652680385221984
    # parameters['spin_period']=10
    # parameters['Wdisk']=1.4460263992995446

    # parameters['dissipation']=True
    # parameters['disk_dissipation_age']=5e-3
    # parameters['wind']=True
    # parameters['wind_saturation_frequency']=2.54
    # parameters['diff_rot_coupling_timescale']=5e-3
    # parameters['wind_strength']=0.17
    # parameters['print_cfile']=False
    # parameters['evolution_max_time_step']=1e-3
    # parameters['evolution_precision']=1e-6
    # parameters['inclination']=0.0
    # parameters['phase_lag_max']= 5.720283932423398e-07
    # parameters['spin_frequency_breaks']=None
    # parameters['spin_frequency_powers']=numpy.array([0.0])
    # parameters['tidal_frequency_breaks']=numpy.array( [0.12566371,0.20936707])
    # parameters['tidal_frequency_powers']=numpy.array( [0.0,0.12528201,0.0])

    #Extremely Slow where Final eccentricity<1e-8
    # parameters['primary_mass']=0.8702751082448668
    # parameters['secondary_mass']=0.8404587704560754
    # parameters['age']=9.168832649763575
    # parameters['feh']= -0.5493213649220865
    # parameters['orbital_period']= 3.42059777531

    # parameters['eccentricity']=0.013787175999258086
    # parameters['spin_period']=10
    # parameters['Wdisk']= 3.9415297170149084

    # parameters['dissipation']=True
    # parameters['disk_dissipation_age']=5e-3
    # parameters['wind']=True
    # parameters['wind_saturation_frequency']=2.54
    # parameters['diff_rot_coupling_timescale']=5e-3
    # parameters['wind_strength']=0.17
    # parameters['print_cfile']=False
    # parameters['evolution_max_time_step']=1e-3
    # parameters['evolution_precision']=1e-6
    # parameters['inclination']=0.0
    # parameters['phase_lag_max']=1.2398415749904416e-06
    # parameters['spin_frequency_breaks']=None
    # parameters['spin_frequency_powers']=numpy.array([0.0])
    # parameters['tidal_frequency_breaks']=numpy.array([0.60887346])
    # parameters['tidal_frequency_powers']=numpy.array([ 0.0, -0.35639464])


    print(parameters)

    evolution=Evolution(interpolator,parameters)
    spin =  evolution.calculate_intial_conditions()
    print('Spin Period = ',spin)







    
    # with open('/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/version2_emcee/catalog/filtering/nominal_value_catalog_temp_cutoff.txt','w') as fnew:
    #     fnew.write('Number\tKIC\tPorb\tspin\teccentricity\tfeh\tage(Gyr)\tm1\tm2\n')
    #     with open('/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/version2_emcee/catalog/filtering/nominal_value_catalog.txt','r') as f:
    #         next(f)
    #         for lines in f:
    #             x=lines.split()
    #             feh=float(x[5])
    #             age=float(x[6])
    #             m1=float(x[7])
    #             m2=float(x[8])
    #             quantity_radius=interpolator('radius',m1, feh)
    #             quantity_lum=interpolator('lum',m1, feh)
    #             t_age=numpy.linspace(5e-3,age,1000)
    #             try:
    #                 T=TeffK(quantity_radius,quantity_lum)(age)
    #                 I=interpolator('iconv',m1,feh)
    #                 if min(I(t_age)):
    #                     fnew.write(lines)
    #                     print(f'System {x[0]} Temp = {T} M1 = {m1} feh = {feh} age = {age}')
    #             except:
    #                 continue


    # I=interpolator('iconv',1.0368188571870431,0.0)

    # t_age=numpy.linspace(5e-3,6,1000)
    # pyplot.plot(t_age,T(t_age))
    # pyplot.show()
