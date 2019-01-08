#!/usr/bin/env python3


import sys
sys.path.append('/Users/ruskinpatel/Desktop/Research/poet/PythonPackage')
sys.path.append('/Users/ruskinpatel/Desktop/Research/poet/scripts')

# import sys
# sys.path.append('.../poet/PythonPackage')
# sys.path.append('.../poet/scripts')


#import sys
#sys.path.append('/home/kpenev/projects/git/poet/PythonPackage')
#sys.path.append('/home/kpenev/projects/git/poet/scripts')


from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from orbital_evolution.binary import Binary
from orbital_evolution.transformations import phase_lag
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet
from mass_calculations import DerivePrimnaryMass
from mass_calculations import DeriveSecondaryMass
from inital_condition_solver import  InitialConditionSolver
from basic_utils import Structure
import numpy
from astropy import units, constants



class evolution:

    def create_planet(self,mass=(constants.M_jup / constants.M_sun).to('')):
        """Return a configured planet to use in the evolution."""

        planet = LockedPlanet(mass=mass, radius=(constants.R_jup / constants.R_sun).to(''))
        return planet

    def create_star(self, mass, dissipation):
        star = EvolvingStar(mass=mass,
                            metallicity=0.0,
                            wind_strength=self.wind_strength if self.wind else 0.0,
                            wind_saturation_frequency=self.wind_saturation_frequency,
                            diff_rot_coupling_timescale=self.diff_rot_coupling_timescale,
                            interpolator=self.interpolator)

        if dissipation == 1:
            star.set_dissipation(zone_index=0,
                                tidal_frequency_breaks=None,
                                spin_frequency_breaks=None,
                                tidal_frequency_powers=numpy.array([0.0]),
                                spin_frequency_powers=numpy.array([0.0]),
                                reference_phase_lag=self.convective_phase_lag)

        #star.set_dissipation(zone_index=1,
        #                     tidal_frequency_breaks=None,
        #                     spin_frequency_breaks=None,
        #                     tidal_frequency_powers=numpy.array([0.0]),
        #                     spin_frequency_powers=numpy.array([0.0]),
        #                     reference_phase_lag=self.)
        return star

    def create_binary_system(self,
                             primary,
                             secondary,
                             initial_semimajor,
                             disk_dissipation_age,
                             secondary_angmom=None):
        """Create a binary system to evolve from the given objects."""

        if isinstance(secondary, LockedPlanet):
            secondary_config = dict(spin_angmom=numpy.array([0.0]),
                                    inclination=None,
                                    periapsis=None)
        else:
           #secondary.select_interpolation_region(disk_dissipation_age)
            secondary_config = dict(spin_angmom=secondary_angmom,
                                    inclination=numpy.array([0.0]),
                                    periapsis=numpy.array([0.0]))

        secondary.configure(age=self.disk_dissipation_age,
                            companion_mass=primary.mass,
                            semimajor=initial_semimajor,
                            eccentricity=0.0,
                            locked_surface=False,
                            zero_outer_inclination=True,
                            zero_outer_periapsis=True,
                            **secondary_config)

        #if isinstance(secondary, EvolvingStar):
            #secondary.detect_stellar_wind_saturation()

        primary.select_interpolation_region(primary.core_formation_age())
        primary.detect_stellar_wind_saturation()

        binary = Binary(primary=primary,
                        secondary=secondary,
                        initial_orbital_period=5.0,
                        initial_eccentricity=0.0,
                        initial_inclination=0.0,
                        disk_lock_frequency=self.disk_lock_frequency,
                        disk_dissipation_age=self.disk_dissipation_age,
                        secondary_formation_age=self.disk_dissipation_age)

        binary.configure(age=primary.core_formation_age(),
                         semimajor=float('nan'),
                         eccentricity=float('nan'),
                         spin_angmom=numpy.array([0.0]),
                         inclination=None,
                         periapsis=None,
                         evolution_mode='LOCKED_SURFACE_SPIN')

        return binary


    def calculate_star_masses(self):

        star_masses = []

        print ("Calculating Masses\n")

        mass1 = DerivePrimnaryMass(
                                    self.interpolator,
                                    self.feh,
                                    self.age,
                                    self.teff)
        PrimaryMass = mass1()


        mass2 = DeriveSecondaryMass(
                                    self.Porb,
                                    self.semimajor,
                                    PrimaryMass)

        SecondaryMass = mass2()

        print ("Primary_Mass = ", PrimaryMass)
        print ("Seconday_Mass = ", SecondaryMass)

        star_masses.append(PrimaryMass)
        star_masses.append(SecondaryMass)

        return star_masses



    def __init__(self,interpolator,observational_parameters,fixed_parameters):

        self.interpolator = interpolator

        self.age = observational_parameters['age']
        self.feh = observational_parameters['feh']
        self.convective_phase_lag = phase_lag(observational_parameters['logQ'])
        self.teff = observational_parameters['teff']
        self.Porb = observational_parameters['Porb']
        self.semimajor = observational_parameters['semimajor']
        self.disk_lock_frequency = observational_parameters['Pdisk']

        self.inclination = fixed_parameters['inclination']
        self.disk_dissipation_age = fixed_parameters['disk_dissipation_age']
        self.wind = fixed_parameters['wind']
        self.planet_formation_age = fixed_parameters['planet_formation_age']
        self.wind_saturation_frequency = fixed_parameters['wind_saturation_frequency']
        self.diff_rot_coupling_timescale = fixed_parameters['diff_rot_coupling_timescale']
        self.wind_strength = fixed_parameters['wind_strength']



    def __call__(self):

        tdisk = self.disk_dissipation_age

        star_masses = self.calculate_star_masses()

        PrimaryMass = star_masses[0]
        SecondaryMass = star_masses[1]

        star = self.create_star(SecondaryMass,0)
        planet = self.create_planet(1.0)

        binary = self.create_binary_system(star,
                                      planet,
                                      10.0,
                                      tdisk)

        binary.evolve(tdisk, 1e-3, 1e-6, None)

        disk_state = binary.final_state()


        planet.delete()
        star.delete()
        binary.delete()

        print ('star-planet evolution completed')

        primary = self.create_star(PrimaryMass, 1)
        secondary = self.create_star(SecondaryMass, 0)
        #secondary = self.create_planet(SecondaryMass)
        find_ic = InitialConditionSolver(disk_dissipation_age=tdisk,
                                         evolution_max_time_step=1.0,
                                         secondary_angmom=numpy.array(
                                             [disk_state.envelope_angmom, disk_state.core_angmom]),
                                         is_secondary_star=True)

        print ('Target parameters: ')
        print ('age = ', self.age)
        print ('Porb = ', self.Porb)
        print ('convective phase lag = ', self.convective_phase_lag)


        target = Structure(age=self.age,
                           Porb=self.Porb,  # current Porb to match
                           Wdisk=self.disk_lock_frequency,  # initial disk locking frequency
                           planet_formation_age=self.planet_formation_age)

        initial_porb, final_Psurf = find_ic(target=target,
                                              primary=primary,
                                              secondary=secondary)

        print ("\nTEST5")


        primary.delete()
        secondary.delete()

        print ("\nTEST6")


        return final_Psurf


