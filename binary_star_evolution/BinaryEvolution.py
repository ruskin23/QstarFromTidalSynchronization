#!/usr/bin/env python3
import shelve
import matplotlib
from matplotlib import pyplot
matplotlib.use('TkAgg')

import sys
sys.path.append('home/ruskin/projects/poet/PythonPackage')
sys.path.append('home/ruskin/projects/poet/scripts')

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

#wsun = 0.24795522138
wsun=1.0

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
                       wsat):
        """plot the evolution of a properly constructed binary."""

        PlotKey=str(self.parameters[self.plot_key])

        if self.plot_primary_envelope==True:pyplot.semilogx(self.Frequencies['age'],
                                                            self.Frequencies['wenv_primary'],
                                                            color='r',#self.plot_color,
                                                            label='Primary Envelope')#_'+self.plot_key+'_'+PlotKey)

        if self.plot_secondary_envelope==True:pyplot.semilogx(self.Frequencies['age'],
                                                            self.Frequencies['wenv_secondary'],
                                                            color='b',#self.plot_color,
                                                            #linestyle=':',
                                                            label='Secondary Envelope')#+self.plot_key+'_'+PlotKey)

        if self.plot_primary_core==True:pyplot.semilogx(self.Frequencies['age'],
                                                        self.Frequencies['wcore_primary'],
                                                        color="r",
                                                        linestyle='--',
                                                        label='Primary Core')#+self.plot_key+'_'+PlotKey)

        if self.plot_secondary_core==True:pyplot.semilogx(self.Frequencies['age'],
                                                          self.Frequencies['wcore_secondary'],
                                                          color="b",
                                                          linestyle='--',
                                                          label='Secondary Core')#+self.plot_key+'_'+PlotKey)



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


        with shelve.open('InitialConditionsFrequencies.shelve','r') as F:
            for key,value in F.items():
                self.Frequencies[key]=value

        if self.plot==True:
            self.plot_evolution(Frequencies,2.54)

    def evolve_binary(self,
                      primary,
                      secondary,
                      SecondaryAngmom):

        binary=self.create_binary_system(primary,
                                         secondary,
                                         secondary_angmom=SecondaryAngmom(self.Wdisk))

        if self.print_cfile==True:
            if self.breaks==True:
                create_c_code='debug/cfile_'+self.system+'_withbreaks.cpp'
            else:
                create_c_code='debug/cfile_'+self.system+'.cpp'

            binary.evolve(
                self.age,
                self.evolution_max_time_step,
                self.evolution_precision,
                None,
                create_c_code=create_c_code,
                eccentricity_expansion_fname=b"eccentricity_expansion_coef.txt")
        else:
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

        self.FinalValues['Porb']=PorbFinal
        self.FinalValues['Eccentricity']=EccFinal
        self.FinalValues['Spin']=spin

        self.Frequencies['age']=evolution.age
        self.Frequencies['wenv_secondary'] = (evolution.secondary_envelope_angmom / binary.secondary.envelope_inertia(evolution.age)) / wsun
        self.Frequencies['wcore_secondary'] = (evolution.secondary_core_angmom / binary.secondary.core_inertia(evolution.age)) / wsun
        self.Frequencies['wenv_primary'] = (evolution.primary_envelope_angmom / binary.primary.envelope_inertia(evolution.age)) / wsun
        self.Frequencies['wcore_primary'] = (evolution.primary_core_angmom / binary.primary.core_inertia(evolution.age)) / wsun
        self.Frequencies['orbitalfrequncy'] = binary.orbital_frequency(evolution.semimajor) / wsun

        tidal_frequency=[]
        for orbf,spinf in zip(self.Frequencies['orbitalfrequncy'],self.Frequencies['wenv_primary']):
            tidal_frequency.append(abs(2*(orbf-spinf)))
        self.Frequencies['tidal_frequency']=tidal_frequency



        if  self.plot==True:

            #pyplot.semilogx(Frequencies['age'],
            #                Frequencies['tidal_frequency'],
            #                color=self.plot_color,
            #                linestyle='-.',
            #                label='tidal_frequency'+self.plot_key+'_'+str(self.parameters[self.plot_key]))

            self.plot_evolution(2.54)

            pyplot.semilogx(self.Frequencies['age'],
                          self.Frequencies['orbitalfrequncy'],
                          linestyle='--',
                          color='k',
                            label='Orbital Frequency')#+self.plot_key+'_'+str(self.parameters[self.plot_key]))



    def __call__(self):

        SecondaryAngmom=IntialSecondaryAngmom(self.interpolator,
                                              self.parameters)

        primary=self.create_star(self.primary_mass)
        secondary=self.create_star(self.secondary_mass)

        if self.GetEvolution==True:
            self.evolve_binary(primary,
                               secondary,
                               SecondaryAngmom)


        if self.GetInitialCondtion==True:
            self.calculate_intial_conditions(primary,
                                             secondary,
                                             SecondaryAngmom)

        if self.ReturnResutls==True:
                return self.FinalValues,self.Frequencies

    def __init__(self,
                 interpolator,
                 parameters):

        self.interpolator=interpolator
        self.parameters=parameters
        for item,value in parameters.items():
            setattr(self,item,value)

        self.FinalValues=dict()
        self.Frequencies=dict()
        self.convective_phase_lag=phase_lag(self.logQ)
        print('Convective Phase lag = ',self.convective_phase_lag)



