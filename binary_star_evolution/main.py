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

wsun = 0.24795522138

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('system',help='pick system from a file')
    parser.add_argument('-b',action='store_const',dest='breaks',
                        const='breaks',
                        help='decide if breaks needed or not')
    args = parser.parse_args()


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
    parameters['PorbInitial']= 9.639020067957873
    parameters['EccentricityInitial']=0.020270562376710374
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

    parameters['GetEvolution']=False
    parameters['GetInitialCondtion']=True


    parameters['plot']=False
    parameters['plot_primary_envelope']=True
    parameters['plot_primary_core']=False
    parameters['plot_secondary_envelope']=False
    parameters['plot_secondary_core']=False
    parameters['plot_key']='logQ'
    parameters['plot_color']='b'

    parameters['logQ']=8.119487908961593

    for key,value in parameters.items():
        print("{} = {}".format(key,value))


    evolve=Evolution(interpolator,parameters)
    evolve()



    #pyplot.axhline(y=(2*numpy.pi/parameters['PspinCurrent'])/wsun,linestyle=':',label='PSpinCurrent',color='r')
    #pyplot.axhline(y=2*numpy.pi/parameters['PorbCurrent']/wsun,linestyle=':',label='PorbCurrent',color='g')
    #pyplot.legend()
    #pyplot.show()

##################################################################################################################################################################################################################################################################################################################################################

    """
    wdisk=[0.8,4.0]
    plot_color=['r','g']
    i=0
    for w in wdisk:
        parameters['Wdisk']=w
        parameters['plot_color']=plot_color[i]
        print('Calculating for {}'.format(w))
        evolve=Evolution(interpolator,parameters)
        evolve()
        i=i+1

    if parameters['plot']==True:

        with shelve.open('InitialConditionsFrequencies.shelve','r') as F:
            age=F['age']
            orbitalfrequncy=F['orbitalfrequncy']

        pyplot.axhline(y=(2*numpy.pi/parameters['PspinCurrent'])/wsun,label='PSpinCurrent',color='r')
        pyplot.axhline(y=2*numpy.pi/parameters['PorbCurrent']/wsun,label='PorbCurrent',color='g')
        pyplot.semilogx(age,orbitalfrequncy,"k",label='OrbitalFrrequency')
        pyplot.legend(loc='upper right')

    pyplot.show()
    """
##################################################################################################################################################################################################################################################################################################################################################





    """
    PorbInitial=[8.3108749748239,8.44329659042446,8.458286134340986,8.459827710559024,8.459982761652514]
    EccentricityInitial=[0.04203515118102658,0.04189305205288432,0.04198917406815888,0.04199891630622043,0.04199989161922022]
    logQ=[6.0,7.0,8.0,9.0,10.0]
    plot_color=["r","g","b","c","m"]

    for i in range(5):

        q=logQ[i]
        print('Caclulating for logQ = ',q)
        parameters['PorbInitial']=PorbInitial[i]
        parameters['EccentricityInitial']=EccentricityInitial[i]
        parameters['plot_color']=plot_color[i]

        evolve=Evolution(interpolator,parameters)
        evolve(q)

    pyplot.show()
    """

