#!/usr/bin/env python3

import matplotlib

matplotlib.use('TkAgg')

import sys

sys.path.append('/Users/ruskinpatel/Desktop/Research/poet/PythonPackage')
sys.path.append('/Users/ruskinpatel/Desktop/Research/poet/scripts')
sys.path.append('/Users/ruskinpatel/Desktop/Research/QstarFromTidalSynchronization')


from matplotlib import pyplot
from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as \
    orbital_evolution_library
from orbital_evolution.binary import Binary
from orbital_evolution.transformations import phase_lag
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet
from basic_utils import Structure
from astropy import units, constants
from basic_utils import Structure
from math import pi
from scipy.optimize import brentq
import numpy

class InitialConditionSolver :
    """Find initial conditions which reproduce a given system now."""

    def _try_initial_conditions(self, initial_orbital_period, disk_period) :
        """
        Get present orbital and stellar spin periods for initial conditions.

        Args:
            - initial_orbital_period:
                The initial orbital period to calculate the deviation for.
            - disk_period:
                The disk locking period to calculate the deviation for.

        Returns:
            - orbital_period:
                The present day orbital period of the system resulting when
                the evolution is started with the input periods.
            - spin_period:
                The present day surface spin of the star resulting when the
                evolution is started with the input periods.
        """

        print('Trying P0 = %s, Pdisk = %s'
              %
              (repr(initial_orbital_period), repr(disk_period)))

        if hasattr(self, 'binary') : self.binary.delete()

        if self.is_secondary_star:
            self.secondary.select_interpolation_region(self.target.disk_dissipation_age)
            secondary_config = dict(spin_angmom=self.secondary_angmom,
                                    inclination=numpy.array([0.0]),
                                    periapsis=numpy.array([0.0]))
            secondary_formation_age = self.target.disk_dissipation_age


        else:
            secondary_config = dict(spin_angmom=numpy.array([0.0]),
                                    inclination=None,
                                    periapsis=None)
            secondary_formation_age = self.target.planet_formation_age




        self.binary = Binary(
            primary = self.primary,
            secondary = self.secondary,
            initial_orbital_period = initial_orbital_period,
            initial_eccentricity = 0.0,
            initial_inclination = 0.0,
            disk_lock_frequency = 2.0 * numpy.pi / disk_period,
            disk_dissipation_age = self.target.disk_dissipation_age,
            secondary_formation_age = secondary_formation_age
        )

        self.binary.primary.select_interpolation_region(self.primary.core_formation_age())
        self.binary.primary.detect_stellar_wind_saturation()


        self.binary.configure(self.primary.core_formation_age(),
                              float('nan'),
                              float('nan'),
                              numpy.array([0.0]),
                              None,
                              None,
                              'LOCKED_SURFACE_SPIN')

        self.binary.secondary.configure(age=self.disk_dissipation_age,
                            companion_mass=self.primary.mass,
                            semimajor=  self.binary.semimajor(initial_orbital_period),
                            eccentricity=0.0,
                            locked_surface=False,
                            zero_outer_inclination=True,
                            zero_outer_periapsis=True,
                            **secondary_config)

        print ("BINARY CONFIGURATION COMPLETE")


        print ("ANGMOM = " , self.primary.mass)

        self.binary.evolve(
            self.target.age,
            self.evolution_max_time_step,
            self.evolution_precision,
            None
        )

        print ("BINARY EVOLUTION COMPLETE")



        final_state = self.binary.final_state()
        assert(final_state.age == self.target.age)
        orbital_period = self.binary.orbital_period(final_state.semimajor)
        stellar_spin_period = (
            2.0 * pi
            *
            self.binary.primary.envelope_inertia(final_state.age)
            /
            final_state.envelope_angmom
        )
        print('Got Porb = %s, P* = %s'
              %
              (repr(orbital_period), repr(stellar_spin_period)))
        if(numpy.isnan(orbital_period)) : orbital_period = 0.0
        return orbital_period, stellar_spin_period


    def __init__(self,
                 planet_formation_age = None,
                 disk_dissipation_age = None,
                 evolution_max_time_step = 1e-3,
                 evolution_precision = 1e-6,
                 orbital_period_tolerance = 1e-6,
                 spin_tolerance = 1e-6,
                 secondary_angmom=None,
                 is_secondary_star = None) :
        """
        Initialize the object.

        Args:
            - planet_formation_age:
                If not None, the planet is assumed to form at the given age
                (in Gyr). Otherwise, the starting age must be specified each
                time this object is called.

            - disk_dissipation_age:
                The age at which the disk dissipates in Gyrs.

            - evolution_max_time_step:
                The maximum timestep the evolution is allowed to make.

            - evolution_precision:
                The precision to require of the evolution.

        Returns: None.
        """

        if planet_formation_age :
            self.planet_formation_age = planet_formation_age
        if disk_dissipation_age is not None :
            self.disk_dissipation_age = disk_dissipation_age
        self.evolution_max_time_step = evolution_max_time_step
        self.evolution_precision = evolution_precision
        self.orbital_period_tolerance = orbital_period_tolerance
        self.spin_tolerance = spin_tolerance
        self.secondary_angmom = secondary_angmom
        self.is_secondary_star = is_secondary_star


    def __call__(self, target, primary, secondary) :

        self.target = target
        self.primary = primary
        self.secondary = secondary

        self._try_initial_conditions(self.target.Porb,self.target.Pdisk)



def create_planet(mass=(constants.M_jup / constants.M_sun).to('')):
    """Return a configured planet to use in the evolution."""

    planet = LockedPlanet(
        mass=mass,
        radius=(constants.R_jup / constants.R_sun).to('')
    )
    return planet


def create_star(mass, interpolator, convective_phase_lag, wind=True):
    star = EvolvingStar(mass=mass,
                        metallicity=0.0,
                        wind_strength=0.17 if wind else 0.0,
                        wind_saturation_frequency=2.78,
                        diff_rot_coupling_timescale=5.0e-3,
                        interpolator=interpolator)

    star.set_dissipation(zone_index=0,
                         tidal_frequency_breaks=None,
                         spin_frequency_breaks=None,
                         tidal_frequency_powers=numpy.array([0.0]),
                         spin_frequency_powers=numpy.array([0.0]),
                         reference_phase_lag=convective_phase_lag)

    star.set_dissipation(zone_index=1,
                         tidal_frequency_breaks=None,
                         spin_frequency_breaks=None,
                         tidal_frequency_powers=numpy.array([0.0]),
                         spin_frequency_powers=numpy.array([0.0]),
                         reference_phase_lag=0.0)
    return star


def create_binary_system(primary,
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

    secondary.configure(age=disk_dissipation_age,
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
                    initial_semimajor=initial_semimajor,
                    initial_eccentricity=0.0,
                    initial_inclination=0.0,
                    disk_lock_frequency=disk_lock_frequency,
                    disk_dissipation_age=disk_dissipation_age,
                    secondary_formation_age=disk_dissipation_age)

    binary.configure(age=primary.core_formation_age(),
                     semimajor=float('nan'),
                     eccentricity=float('nan'),
                     spin_angmom=numpy.array([0.0]),
                     inclination=None,
                     periapsis=None,
                     evolution_mode='LOCKED_SURFACE_SPIN')

    return binary



def test_ic_solver(interpolator,convective_phase_lag,wind):
    """Find initial condition to reproduce some current state and plot."""

    tdisk = 5e-3


    star = create_star(0.8, interpolator=interpolator, convective_phase_lag=0.0, wind=wind)
    planet = create_planet(1.0)

    binary = create_binary_system(star,
                                  planet,
                                  2.0 * numpy.pi / 3.0,
                                  10.0,
                                  tdisk)

    binary.evolve(tdisk, 1e-3, 1e-6, None)

    disk_state = binary.final_state()

    print("FINISHED PLANET-STAR EVOLUTION")

    planet.delete()
    star.delete()
    binary.delete()

    find_ic = InitialConditionSolver(disk_dissipation_age=5e-3,
                                    evolution_max_time_step=1e-2,
                                    secondary_angmom=numpy.array([disk_state.envelope_angmom, disk_state.core_angmom]),
                                    is_secondary_star=True)

    P_disk = 1.0
    P_spin = []


    primary = create_star(1.0, interpolator, convective_phase_lag, wind = wind)
    secondary = create_star(0.8, interpolator, convective_phase_lag, wind = wind)


    target = Structure( age=7.0,
                        Porb=3.0, #initially 3.0
                        Pdisk=1.0, #initially Psurf=10.0
                        planet_formation_age=5e-3)

    find_ic(target=target, primary=primary, secondary=secondary)



    primary.delete()
    secodnary.delete()
    binary.delete()


if __name__ == '__main__':
    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        b"eccentricity_expansion_coef.txt"
    )
    serialized_dir = '/Users/ruskinpatel/Desktop/Research/poet/stellar_evolution_interpolators'
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    # test_evolution(interpolator, phase_lag(6.0))
    test_ic_solver(interpolator,phase_lag(6.0),True)
