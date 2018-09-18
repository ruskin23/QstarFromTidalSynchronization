#!/usr/bin/env python3

import matplotlib
matplotlib.use('TkAgg')

import sys
sys.path.append('/Users/ruskinpatel/Desktop/Research/poet/PythonPackage')
sys.path.append('/Users/ruskinpatel/Desktop/Research/poet/scripts')

from matplotlib import pyplot
from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from orbital_evolution.binary import Binary
from orbital_evolution.transformations import phase_lag
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet
import numpy
from astropy import units, constants



class binary_evolution:

    def create_planet(self,mass=(constants.M_jup / constants.M_sun).to('')):
        """Return a configured planet to use in the evolution."""

        planet = LockedPlanet(mass=mass, radius=(constants.R_jup / constants.R_sun).to(''))
        return planet

    def create_star(self,mass):
        star = EvolvingStar(mass=mass,
                            metallicity=0.0,
                            wind_strength=0.17 if self.wind else 0.0,
                            wind_saturation_frequency=self.wind_saturation_frequency,
                            diff_rot_coupling_timescale=self.diff_rot_coupling_timescale,
                            interpolator=self.interpolator)

        star.set_dissipation(zone_index=0,
                             tidal_frequency_breaks=None,
                             spin_frequency_breaks=None,
                             tidal_frequency_powers=numpy.array([0.0]),
                             spin_frequency_powers=numpy.array([0.0]),
                             reference_phase_lag=self.convective_phase_lag)

        star.set_dissipation(zone_index=1,
                             tidal_frequency_breaks=None,
                             spin_frequency_breaks=None,
                             tidal_frequency_powers=numpy.array([0.0]),
                             spin_frequency_powers=numpy.array([0.0]),
                             reference_phase_lag=0.0)
        return star

    def create_binary_system(self,
                             primary,
                             secondary,
                             disk_lock_frequency,
                             initial_semimajor,
                             disk_dissipation_age,
                             secondary_angmom=None):
        """Create a binary system to evolve from the given objects."""

        if isinstance(secondary, LockedPlanet):
            secondary_config = dict(spin_angmom=numpy.array([0.0]),
                                    inclination=None,
                                    periapsis=None)
        else:
            secondary.select_interpolation_region(disk_dissipation_age)
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

        if isinstance(secondary, EvolvingStar):
            secondary.detect_stellar_wind_saturation()

        primary.select_interpolation_region(primary.core_formation_age())
        primary.detect_stellar_wind_saturation()

        binary = Binary(primary=primary,
                        secondary=secondary,
                        initial_orbital_period=initial_orbital_period,
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

        mass1 = DerivePrimnaryMass(self.interpolator, self.feh, self.age, self.teff)
        mass2 = DeriveSecondaryMass(self.Porb,self.rvk,self.inclination)

        PrimaryMass = mass1()
        SecondaryMass = mass2(PrimaryMass)

        star_masses.append(PrimaryMass)
        star_masses.append(SecondaryMass)

        return star_masses



    def __init__(self,interpolator,parameters,fixed_parameters):

        self.interpolator = interpolator
        self.age = parameters['age']
        self.convective_phase_lag = parameters['convective_phase_lag']
        self.teff = parameters['teff']
        self.Porb = parameters['Porb']
        self.rvk = parameters['rvk']
        self.inlination = fixed_parameters['inclination']
        self.disk_dissipation_age = fixed_parameters['disk_dissipation_age']
        self.wind = fixed_parameters['wind']
        self.planet_formation_age = fixed_parameters['planet_formation_age']
        self.disk_lock_frequency = fixed_parameters['disk_lock_frequency']
        self.wind_saturation_frequency = fixed_parameters['wind_saturation_frequency']
        self.diff_rot_coupling_timescale = fixed_parameters['diff_rot_coupling_timescale']



    def __call__(self):

        tdisk = self.disk_dissipation_age

        star_masses = self.calculate_star_masses()

        PrimaryMass = star_masses[0]
        SecondaryMass = star_masses[1]

        star = self.create_star(SecondaryMass)
        planet = self.create_planet(1.0)

        binary = create_binary_system(star,
                                      planet,
                                      2.0 * numpy.pi / 3.0,
                                      10.0,
                                      tdisk)

        binary.evolve(tdisk, 1e-3, 1e-6, None)

        disk_state = binary.final_state()

        planet.delete()
        star.delete()
        binary.delete()

        primary = self.create_star(PrimaryMass)
        secondary = self.create_star(SecondaryMass)

        find_ic = InitialConditionSolver(disk_dissipation_age=tdisk,
                                         evolution_max_time_step=1e-2,
                                         secondary_angmom=numpy.array(
                                             [disk_state.envelope_angmom, disk_state.core_angmom]),
                                         is_secondary_star=True)

        target = Structure(age=self.age,
                           Porb=self.Porb,  # current Porb to match
                           Wdisk=self.disk_lock_frequency,  # initial disk locking frequency
                           planet_formation_age=self.planet_formation_age)

        initial_porb, final_Psurf = find_ic(target=target,
                                              primary=primary,
                                              secondary=secondary)

        primary.delete()
        secondary.delete()
        binary.delete()

        return final_Psurf


