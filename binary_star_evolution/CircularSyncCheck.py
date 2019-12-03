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

    system=args.system

    with open('SpinlogQCatalog_el0.4.txt','r') as f:
        for lines in f:
            x=lines.split()
            if x[0]==system:
                primary_mass=float(x[15])
                massratio=float(x[14])
                secondary_mass=massratio*primary_mass
                age=float(x[16])
                feh=float(x[17])
                PspinCurrent=float(x[12])
                PorbCurrent=float(x[6])
                EccentricityCurrent=float(x[8])

    parameters=dict()

    parameters['system']=system
    parameters['primary_mass']=primary_mass
    parameters['secondary_mass']= secondary_mass
    parameters['age']=age
    parameters['feh']=feh
    parameters['PorbInitial']= PorbCurrent
    parameters['EccentricityInitial']=EccentricityCurrent
    parameters['PorbCurrent']=PorbCurrent
    parameters['PspinCurrent']=PspinCurrent
    parameters['EccentricityCurrent']=EccentricityCurrent

    if args.breaks:
        #TidalFrequencyBreaks=numpy.array([abs(-4*numpy.pi*((1.0/parameters['PorbCurrent'])
        #                                                   -
        #                                                   1.0/parameters['PspinCurrent']))])

        TidalFrequencyBreaks=numpy.array([abs(2*numpy.pi)])
        TidalFrequencyPowers=numpy.array([2.0,2.0])
        parameters['breaks']=True
    else:
        TidalFrequencyBreaks=None
        TidalFrequencyPowers=numpy.array([0.0])
        parameters['breaks']=False

    parameters['tidal_frequency_breaks']=TidalFrequencyBreaks
    parameters['tidal_frequency_powers']=TidalFrequencyPowers

    parameters['spin_frequency_breaks']=None
    parameters['spin_frequency_powers']=numpy.array([0.0])

    parameters['dissipation']=True

    parameters['Wdisk']=4.1
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

    parameters['plot']=True
    parameters['plot_primary_envelope']=True
    parameters['plot_primary_core']=False
    parameters['plot_secondary_envelope']=False
    parameters['plot_secondary_core']=False
    parameters['plot_key']='logQ'
    parameters['plot_color']='b'

    parameters['logQ']=6.0

    for key,value in parameters.items():
        print("{} = {}".format(key,value))

