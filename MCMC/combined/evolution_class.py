#!/usr/bin/env python3 -u

import sys

from pathlib import Path
home_dir=str(Path.home())

if home_dir=='/home/rxp163130':
    poet_path=home_dir+'/poet'
if home_dir=='/home/ruskin':
    poet_path=home_dir+'/projects/poet'

sys.path.append(poet_path+'/PythonPackage')
sys.path.append(poet_path+'/scripts')

import sys

from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from orbital_evolution.binary import Binary
from orbital_evolution.transformations import phase_lag
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet
from stellar_evolution.manager import StellarEvolutionManager
from stellar_evolution.derived_stellar_quantities import\
    TeffK,\
    LogGCGS,\
    RhoCGS
from initial_conditions_solver import  InitialConditionSolver
from basic_utils import Structure
import numpy
import scipy
from astropy import units, constants
import pickle


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




    def __init__(self,
                 interpolator,
                 sampling_parameters,
                 fixed_parameters,
                 mass_ratio,
                 instance,
                 output_directory):

        self.interpolator=interpolator

        for item,value in sampling_parameters.items():
            setattr(self,item,value)

        for item,value in fixed_parameters.items():
            setattr(self,item,value)


        self.convective_phase_lag=phase_lag(self.logQ)

        self.mass_ratio=mass_ratio
        self.instance=instance
        self.output_directory=output_directory

        self.model_parameters=dict()

    def __call__(self):

        tdisk = self.disk_dissipation_age

        self.secondary_mass=self.primary_mass*self.mass_ratio
        if numpy.logical_or((numpy.logical_or(self.secondary_mass>1.2,
                                              self.secondary_mass<0.4)),
                            (numpy.logical_or(numpy.isnan(self.primary_mass),
                                              numpy.isnan(self.secondary_mass)))
                            ):
            print('mass out of range')
            sys.stdout.flush()
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
        sys.stdout.flush()
        primary = self.create_star(self.primary_mass, 1)
        secondary = self.create_star(self.secondary_mass, 1)
        find_ic = InitialConditionSolver(disk_dissipation_age=tdisk,
                                         evolution_max_time_step=1e-3,
                                         secondary_angmom=numpy.array(
                                             [disk_state.envelope_angmom, disk_state.core_angmom]),
                                         is_secondary_star=True,
                                         instance = self.instance,
                                         output_directory=self.output_directory)

        target = Structure(age=self.age,
                           Porb=self.Porb,  # current Porb to match
                           Wdisk=self.Wdisk,
                           eccentricity=self.eccentricity)

        self.spin=find_ic(target=target,
                          primary=primary,
                          secondary=secondary)



        primary.delete()
        secondary.delete()


        return self.spin
