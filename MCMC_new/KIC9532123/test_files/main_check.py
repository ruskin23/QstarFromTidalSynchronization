#!/usr/bin/env python3


import sys
sys.path.append('.../poet/PythonPackage')
sys.path.append('.../poet/scripts')

from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from orbital_evolution.binary import Binary
from orbital_evolution.transformations import phase_lag
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet
from mass_calculation import Derive_mass
from solver_check import  InitialConditionSolver
from basic_utils import Structure
import numpy
import scipy
from astropy import units, constants
import pickle
import argparse



class evolution:

    def create_planet(self,mass=(constants.M_jup / constants.M_sun).to('')):
        """Return a configured planet to use in the evolution."""

        planet = LockedPlanet(mass=mass, radius=(constants.R_jup / constants.R_sun).to(''))
        return planet

    def create_star(self, mass, dissipation):
        star = EvolvingStar(mass=mass,
                            metallicity=self.feh,
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

        print ("BEGINSAT")

        if isinstance(secondary, EvolvingStar):
            secondary.detect_stellar_wind_saturation()
            print ("DETECTED SECONDARY WIND SAT")



        primary.select_interpolation_region(primary.core_formation_age())
        primary.detect_stellar_wind_saturation()

        binary = Binary(primary=primary,
                        secondary=secondary,
                        initial_orbital_period=self.Porb,
                        initial_eccentricity=0.0,
                        initial_inclination=0.0,
                        disk_lock_frequency=self.Wdisk,
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

        mass_ratio = 0.8840  # mass ratio is M2/M1

        mass = Derive_mass(
                            self.interpolator,
                            self.teff_primary,
                            self.logg,
                            self.feh)

        sol = mass()
        self.primary_mass = sol[0]
        self.age=sol[1]
        self.secondary_mass = self.primary_mass*mass_ratio

        print(self.primary_mass)
        print(self.secondary_mass)
        print(self.age)


    def __init__(self,interpolator,observational_parameters):

        self.interpolator=interpolator

        for item,value in observational_parameters.items():
            setattr(self,item,value)

        self.convective_phase_lag=phase_lag(self.logQ)

        self.primary_mass=0.0
        self.secondary_mass=0.0



    def __call__(self):

        tdisk = self.disk_dissipation_age

        self.calculate_star_masses()
        if numpy.logical_or( numpy.isnan(self.primary_mass),
                            numpy.isnan(self.secondary_mass)):
            print('mass out of range')
            return scipy.nan

        star = self.create_star(self.secondary_mass,1)
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

        primary = self.create_star(self.primary_mass, 1)
        secondary = self.create_star(self.secondary_mass, 1)
        find_ic = InitialConditionSolver(primary,
                                         secondary,
                                         disk_dissipation_age=tdisk,
                                         evolution_max_time_step=1e-3,
                                         secondary_angmom=numpy.array([disk_state.envelope_angmom, disk_state.core_angmom]),
                                         is_secondary_star=True)

        target = Structure(age=self.age,
                           Porb=self.Porb,  # current Porb to match
                           Wdisk=self.Wdisk,
                           eccentricity=self.eccentricity,
                           logQ=self.logQ)

        find_ic(target=target)

        primary.delete()
        secondary.delete()


if __name__ == '__main__':


    for i in range(1):
        serialized_dir ="/home/ruskin/projects/poet/stellar_evolution_interpolators"
        manager = StellarEvolutionManager(serialized_dir)
        interpolator = manager.get_interpolator_by_name('default')

        orbital_evolution_library.read_eccentricity_expansion_coefficients(
            b"eccentricity_expansion_coef.txt"
        )

        parameters = dict(
                      teff_primary= 5335.492032128414,
                      feh=-0.7690531666974676,
                      logg=4.621887887981691,
                      logQ=5.332,
                      Wdisk=4.298549388598006,
                      Porb=7.903632794647903,
                      incination=0.0,
                      disk_dissipation_age=5e-3,
                      wind=True,
                      planet_formation_age=5e-3,
                      wind_saturation_frequency=2.54,
                      diff_rot_coupling_timescale=5e-3,
                      wind_strength=0.17,
                      eccentricity=0.20782473351422318)


        evolve_check = evolution(interpolator,parameters)
        evolve_check()




