import sys
sys.path.append('home/ruskin/projects/poet/PythonPackage')
sys.path.append('home/ruskin/projects/poet/scripts')

from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from orbital_evolution.binary import Binary
from orbital_evolution.transformations import phase_lag
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet
from basic_utils import Structure
import numpy
import scipy
from astropy import units, constants


class IntialSecondaryAngmom:

    def create_planet(self,mass=(constants.M_jup / constants.M_sun).to('')):
        """Return a configured planet to use in the evolution."""

        planet = LockedPlanet(mass=mass, radius=(constants.R_jup / constants.R_sun).to(''))
        return planet

    def create_star(self, mass):
        star = EvolvingStar(mass=mass,
                            metallicity=self.feh,
                            wind_strength=self.wind_strength if self.wind else 0.0,
                            wind_saturation_frequency=self.wind_saturation_frequency,
                            diff_rot_coupling_timescale=self.diff_rot_coupling_timescale,
                            interpolator=self.interpolator)

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





    def __init__(self,
                 interpolator,
                 logQ,
                 parameters
                 ):


        self.interpolator=interpolator
        self.convective_phase_lag=phase_lag(logQ)

        for item,value in parameters.items():
            setattr(self,item,value)


    def __call__(self,wdisk):


        self.disk_lock_frequency=wdisk
        star = self.create_star(self.secondary_mass)
        planet = self.create_planet(1.0)

        binary = self.create_binary_system(star,
                                      planet,
                                      10.0,
                                      self.disk_dissipation_age)

        binary.evolve(self.disk_dissipation_age, 1e-3, 1e-6, None)
        disk_state = binary.final_state()

        planet.delete()
        star.delete()
        binary.delete()

        print ('star-planet evolution completed')

        return numpy.array([disk_state.envelope_angmom, disk_state.core_angmom])

