#!/usr/bin/env python3 -u

import sys

sys.path.append('/home/ruskin/projects/poet/PythonPackage')
sys.path.append('/home/ruskin/projects/poet/scripts')

from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from orbital_evolution.binary import Binary
from orbital_evolution.transformations import phase_lag
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet

import numpy
import scipy
from astropy import units, constants


class BinaryObjects():

    def create_planet(self,mass=(constants.M_jup / constants.M_sun).to('')): 
        """Return a configured planet to use in the evolution."""

        planet = LockedPlanet(mass=mass, radius=(constants.R_jup / constants.R_sun).to(''))
        return planet

    def create_star(self, mass, dissipation=None):
        star = EvolvingStar(mass=mass,
                            metallicity=self.feh,
                            wind_strength=self.wind_strength if self.wind else 0.0, 
                            wind_saturation_frequency=self.wind_saturation_frequency, 
                            diff_rot_coupling_timescale=self.diff_rot_coupling_timescale, 
                            interpolator=self.interpolator)

        if dissipation is True:
            star.set_dissipation(zone_index=0,
                                tidal_frequency_breaks=self.tidal_frequency_breaks,
                                spin_frequency_breaks=self.spin_frequency_breaks,
                                tidal_frequency_powers=self.tidal_frequency_powers,
                                spin_frequency_powers=self.spin_frequency_powers,
                                reference_phase_lag=self.convective_phase_lag)


        return star


    def create_binary_system(self,
                             primary,
                             secondary,
                             secondary_angmom=None):
        """Create a binary system to evolve from the given objects."""

        if isinstance(secondary, LockedPlanet):
            spin_angmom=numpy.array([0.0])
            inclination=None
            periapsis=None
        else:
            secondary.select_interpolation_region(self.disk_dissipation_age)
            spin_angmom=secondary_angmom
            inclination=numpy.array([0.0])
            periapsis=numpy.array([0.0])

        binary = Binary(primary=primary,
                        secondary=secondary,
                        initial_orbital_period=self.orbital_period,
                        initial_eccentricity=self.eccentricity,
                        initial_inclination=0.0,
                        disk_lock_frequency=self.Wdisk, 
                        disk_dissipation_age=self.disk_dissipation_age, 
                        secondary_formation_age=self.disk_dissipation_age) 
        
        binary.primary.select_interpolation_region(primary.core_formation_age())
        if isinstance(secondary, EvolvingStar):binary.secondary.detect_stellar_wind_saturation()

        binary.configure(age=primary.core_formation_age(),
                         semimajor=float('nan'),
                         eccentricity=float('nan'),
                         spin_angmom=numpy.array([0.0]),
                         inclination=None,
                         periapsis=None,
                         evolution_mode='LOCKED_SURFACE_SPIN')

        binary.primary.detect_stellar_wind_saturation()

        binary.secondary.configure(age=self.disk_dissipation_age, 
                            companion_mass=primary.mass,
                            semimajor=binary.semimajor(self.orbital_period),
                            eccentricity=self.eccentricity,
                            spin_angmom=spin_angmom,
                            inclination=inclination,
                            periapsis=periapsis,
                            locked_surface=False,
                            zero_outer_inclination=True,
                            zero_outer_periapsis=True
                            )


        return binary


    def __init__(self,
                interpolator,
                parameters):

        self.interpolator=interpolator

        for item,value in parameters.items():
            setattr(self,item,value)
        self.convective_phase_lag=phase_lag(self.logQ)
