import numpy
import argparse
import matplotlib
from matplotlib import pyplot
matplotlib.use('TkAgg')

import sys
sys.path.append('home/ruskin/projects/poet/PythonPackage')
sys.path.append('home/ruskin/projects/poet/scripts')

from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as \
    orbital_evolution_library
from BinaryEvolution import Evolution




if __name__=='__main__':

    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        b"eccentricity_expansion_coef.txt"
    )
    serialized_dir = '/home/kpenev/projects/git/poet/stellar_evolution_interpolators'

    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    parameters=dict()

    primary_masses=[0.4,0.6,0.8,1.0,1.2]
    secondary_masses=[1.2,1.0,0.8,0.6,0.4]
    wdisk=[0.5,1.0,2.0,3.0,4.0]

    parameters['age']=6.0
    parameters['feh']=0.0


    parameters['PorbInitial']= 8.0
    parameters['EccentricityInitial']=0.1

    TidalFrequencyBreaks=numpy.array([abs(2*numpy.pi)])
    TidalFrequencyPowers=numpy.array([2.0,2.0])
    parameters['breaks']=True

    parameters['tidal_frequency_breaks']=TidalFrequencyBreaks
    parameters['tidal_frequency_powers']=TidalFrequencyPowers

    parameters['spin_frequency_breaks']=None
    parameters['spin_frequency_powers']=numpy.array([0.0])

    parameters['dissipation']=True

    parameters['disk_dissipation_age']=5e-3
    parameters['wind']=True
    parameters['wind_saturation_frequency']=2.54
    parameters['diff_rot_coupling_timescale']=5e-3
    parameters['wind_strength']=0.17

    parameters['evolution_max_time_step']=1e-3
    parameters['evolution_precision']=1e-6
    parameters['print_cfile']=False

    parameters['GetEvolution']=True
    parameters['GetInitialCondtion']=False

    parameters['ReturnResults']=True

    parameters['plot']=False
    parameters['plot_primary_envelope']=True
    parameters['plot_primary_core']=False
    parameters['plot_secondary_envelope']=False
    parameters['plot_secondary_core']=False
    parameters['plot_key']='logQ'
    parameters['plot_color']='b'


    for key,value in parameters.items():
        print("{} = {}".format(key,value))

    logQ=[3.0,3.5,4.0,4.5,5.0]

    with open('CheckCircularVsSync.txt','w') as f:
        f.write('PrimaryMass'+'\t'+
                'SecondaryMass'+'\t'+
                'Wdisk'+'\t'+
                'logQ'+'\t'+
                'FinalPorb'+'\t'+
                'FinalEccentricity'+'\t'+
                'FinalSpin'+'\n')

    for i in range(5):
        FinalResults=dict()

        parameters['primary_mass']=primary_masses[i]
        parameters['secondary_mass']=secondary_masses[i]
        parameters['Wdisk']=wdisk[i]
        parameters['logQ']=logQ[i]

        print('Calculating for paramters:',parameters)


        evolve=Evolution(interpolator,parameters)
        FinalResults=evolve()


        with open('CheckCircularVsSync.txt','a') as f:
            f.write(repr(parameters['primary_mass'])+'\t'+
                    repr(parameters['secondary_mass'])+'\t'+
                    repr(parameters['Wdisk'])+'\t'+
                    repr(parameters['logQ'])+'\t'+
                    repr(FinalResults['Porb'])+'\t'+
                    repr(FinalResults['Eccentricity'])+'\t'+
                    repr(FinalResults['spin'])+'\n')
