#!/usr/bin/env python3

import argparse
import shelve
import matplotlib

matplotlib.use('TkAgg')

import sys
sys.path.append('home/ruskin/projects/poet/PythonPackage')
sys.path.append('home/ruskin/projects/poet/scripts')


from matplotlib import pyplot
from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as \
    orbital_evolution_library
from orbital_evolution.binary import Binary
from orbital_evolution.transformations import phase_lag
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet
from InitialSecondaryAngmom import IntialSecondaryAngmom
from initial_condition_solver import  InitialConditionSolver
from basic_utils import Structure
import numpy
import scipy
from astropy import units, constants

wsun = 0.24795522138

class Evolution:

    def create_star(self,
                    mass):

        star = EvolvingStar(mass=mass,
                            metallicity=self.feh,
                            wind_strength=self.wind_strength if self.wind else 0.0,
                            wind_saturation_frequency=self.wind_saturation_frequency,
                            diff_rot_coupling_timescale=self.diff_rot_coupling_timescale,
                            interpolator=self.interpolator)

        if self.dissipation == True:
            star.set_dissipation(zone_index=0,
                                 tidal_frequency_breaks=self.tidal_frequency_breaks,
                                 spin_frequency_breaks=self.spin_frequency_breaks,
                                 tidal_frequency_powers=self.tidal_frequency_powers,
                                 spin_frequency_powers=self.spin_frequency_powers,
                                 reference_phase_lag=self.convective_phase_lag)

            #star.set_dissipation(zone_index=1,
            #                     tidal_frequency_breaks=None,
            #                    spin_frequency_breaks=None,
            #                     tidal_frequency_powers=numpy.array([0.0]),
            #                     spin_frequency_powers=numpy.array([0.0]),
            #                     reference_phase_lag=convective_phase_lag)

        return star



    def create_binary_system(self,
                             primary,
                             secondary,
                             secondary_angmom=None):

        if isinstance(secondary, LockedPlanet):
            secondary_config = dict(spin_angmom=numpy.array([0.0]),
                                    inclination=None,
                                    periapsis=None)
        else:
            secondary.select_interpolation_region(self.disk_dissipation_age)
            secondary_config = dict(spin_angmom=secondary_angmom,
                                    inclination=numpy.array([0.0]),
                                    periapsis=numpy.array([0.0]))
        print(secondary_config['spin_angmom'])

        binary = Binary(primary=primary,
                        secondary=secondary,
                        initial_orbital_period=self.PorbInitial,
                        initial_eccentricity=self.EccentricityInitial,
                        initial_inclination=0.0,
                        disk_lock_frequency=self.Wdisk,
                        disk_dissipation_age=self.disk_dissipation_age,
                        secondary_formation_age=self.disk_dissipation_age)


        primary.select_interpolation_region(primary.core_formation_age())
        if isinstance(secondary, EvolvingStar):secondary.detect_stellar_wind_saturation()


        binary.configure(age=primary.core_formation_age(),
                         semimajor=float('nan'),
                         eccentricity=float('nan'),
                         spin_angmom=numpy.array([0.0]),
                         inclination=None,
                         periapsis=None,
                         evolution_mode='LOCKED_SURFACE_SPIN')

        primary.detect_stellar_wind_saturation()

        binary.secondary.configure(age=self.disk_dissipation_age,
                            companion_mass=primary.mass,
                            semimajor=binary.semimajor(self.PorbInitial),
                            eccentricity=self.EccentricityInitial,
                            locked_surface=False,
                            zero_outer_inclination=True,
                            zero_outer_periapsis=True,
                            **secondary_config)

        return binary


    def plot_evolution(self,
                       Frequencies,
                       wsat):
        """plot the evolution of a properly constructed binary."""

        PlotKey=str(parameters[self.plot_key])

        if self.plot_primary_envelope==True:pyplot.semilogx(Frequencies['age'],
                                                            Frequencies['wenv_primary'],
                                                            color=self.plot_color,
                                                            label='PrimaryEnvelope_'+self.plot_key+'_'+PlotKey)

        if self.plot_secondary_envelope==True:pyplot.semilogx(Frequencies['age'],
                                                            Frequencies['wenv_secondary'],
                                                            color=self.plot_color,
                                                            linestyle=':',
                                                            label='SecondaryEnvelope_logQ'+self.plot_key+'_'+PlotKey)

        if self.plot_primary_core==True:pyplot.semilogx(Frequencies['age'],
                                                        Frequencies['wcore_primary'],
                                                        color="b",
                                                        linestyle='--',
                                                        label='PrimaryCore_logQ'+self.plot_key+'_'+PlotKey)

        if self.plot_secondary_core==True:pyplot.semilogx(Frequencies['age'],
                                                          Frequencies['wcore_secondary'],
                                                          color="r",linestyle='--',
                                                          label='SecondaryCore_logQ'+self.plot_key+'_'+PlotKey)



    def calculate_intial_conditions(self,
                              primary,
                              secondary,
                              SecondaryAngmom):

        FindIC=InitialConditionSolver(system=self.system,
                                      print_cfile=self.print_cfile,
                                      breaks=self.breaks,
                                      disk_dissipation_age=self.disk_dissipation_age,
                                      evolution_max_time_step=1e-3,
                                      secondary_angmom=SecondaryAngmom(self.Wdisk),
                                      is_secondary_star=True)

        target=Structure(age=self.age,
                         Porb=self.PorbCurrent,
                         Wdisk=self.Wdisk,
                         eccentricity=self.EccentricityCurrent)

        Solutions=FindIC(target=target,
                         primary=primary,
                         secondary=secondary)

        print('Final Pspin = ', Solutions['spin'])

        Frequencies=dict()

        if self.plot==True:
            with shelve.open('InitialConditionsFrequencies.shelve','r') as F:
                for key,value in F.items():
                    Frequencies[key]=value

            self.plot_evolution(Frequencies,2.54)

    def evolve_binary(self,
                      primary,
                      secondary,
                      SecondaryAngmom):

        binary=self.create_binary_system(primary,
                                         secondary,
                                         secondary_angmom=SecondaryAngmom(self.Wdisk))
        binary.evolve(self.age,
                      self.evolution_max_time_step,
                      self.evolution_precision,
                      None)

        final_state=binary.final_state()
        evolution=binary.get_evolution()

        spin = (2.0*numpy.pi*binary.primary.envelope_inertia(final_state.age)/
                final_state.primary_envelope_angmom)

        PorbFinal=binary.orbital_period(final_state.semimajor)
        EccFinal=final_state.eccentricity

        print('Final Porb = ',PorbFinal )
        print('Final Eccentricity = ', EccFinal)
        print('Secondary Initial Angmom = ', SecondaryAngmom(self.Wdisk))
        print('Final Spin = ',spin)

        if  self.plot==True:

            Frequencies=dict()

            Frequencies['age']=evolution.age
            Frequencies['wenv_secondary'] = (evolution.secondary_envelope_angmom / binary.secondary.envelope_inertia(evolution.age)) / wsun
            Frequencies['wcore_secondary'] = (evolution.secondary_core_angmom / binary.secondary.core_inertia(evolution.age)) / wsun
            Frequencies['wenv_primary'] = (evolution.primary_envelope_angmom / binary.primary.envelope_inertia(evolution.age)) / wsun
            Frequencies['wcore_primary'] = (evolution.primary_core_angmom / binary.primary.core_inertia(evolution.age)) / wsun
            Frequencies['orbitalfrequncy'] = binary.orbital_frequency(evolution.semimajor) / wsun

            self.plot_evolution(Frequencies,2.54)


    def __call__(self):

        self.convective_phase_lag=phase_lag(self.logQ)
        SecondaryAngmom=IntialSecondaryAngmom(self.interpolator,
                                              self.parameters)

        primary=self.create_star(self.primary_mass)
        secondary=self.create_star(self.secondary_mass)

        if self.GetEvolution==True:self.evolve_binary(primary,
                                                      secondary,
                                                      SecondaryAngmom)


        if self.GetInitialCondtion==True:self.calculate_intial_conditions(primary,
                                                                          secondary,
                                                                          SecondaryAngmom)


    def __init__(self,
                 interpolator,
                 parameters):

        self.interpolator=interpolator
        self.parameters=parameters
        for item,value in parameters.items():
            setattr(self,item,value)

if __name__ == '__main__':

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

    with open('spin_vs_logQ_systems_0.2.txt','r') as f:
        for lines in f:
            x=lines.split()
            if x[0]==system:
                primary_mass=float(x[15])
                massratio=float(x[14])
                secondary_mass=massratio*primary_mass
                age=10**(float(x[16]))
                feh=float(x[17])
                PspinCurrent=float(x[12])
                PorbCurrent=float(x[6])
                EccentricityCurrent=float(x[8])

    parameters=dict()

    parameters['system']=system
    parameters['primary_mass']=primary_mass
    parameters['secondary_mass']=secondary_mass
    parameters['age']=age
    parameters['feh']=feh
    parameters['PorbInitial']=PorbCurrent
    parameters['EccentricityInitial']=EccentricityCurrent
    parameters['PorbCurrent']=PorbCurrent
    parameters['PspinCurrent']=PspinCurrent
    parameters['EccentricityCurrent']=EccentricityCurrent

    if args.breaks:
        TidalFrequencyBreaks=numpy.array([abs(-4*numpy.pi*((1.0/parameters['PorbCurrent'])
                                                           -
                                                           1.0/parameters['PspinCurrent']))])
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

    #parameters['Wdisk']=4.1
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


    parameters['plot']=True
    parameters['plot_primary_envelope']=True
    parameters['plot_primary_core']=False
    parameters['plot_secondary_envelope']=False
    parameters['plot_secondary_core']=False
    parameters['plot_key']='Wdisk'
    #parameters['plot_color']='-b'

    parameters['logQ']=6.021021021021021

    for key,value in parameters.items():
        print("{} = {}".format(key,value))

    wdisk=[0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0]
    plot_color=['r','g','b','c','m','y','k','tab:olive']
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

        pyplot.axhline(y=(2*numpy.pi/parameters['PspinCurrent'])/wsun,label='PSpinCurrent')
        pyplot.axhline(y=2*numpy.pi/parameters['PorbCurrent']/wsun,label='PorbCurrent')
        pyplot.semilogx(age,orbitalfrequncy,"k",label='OrbitalFrrequency')
        pyplot.legend(loc='upper right')

    pyplot.show()

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
