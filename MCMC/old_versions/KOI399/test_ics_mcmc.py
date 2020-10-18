import sys
sys.path.append('.../poet/PythonPackage')
sys.path.append('.../poet/scripts')


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
import scipy
import numpy

from inital_condition_solver import  InitialConditionSolver


class ics_check():

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

        print ("BEGINSAT")

        if isinstance(secondary, EvolvingStar):
            secondary.detect_stellar_wind_saturation()
            print ("DETECTED SECONDARY WIND SAT")




        primary.select_interpolation_region(primary.core_formation_age())
        primary.detect_stellar_wind_saturation()

        binary = Binary(primary=primary,
                        secondary=secondary,
                        initial_orbital_period=self.initial_orbital_period,
                        initial_eccentricity=0.0,
                        initial_inclination=0.0,
                        disk_lock_frequency=self.disk_lock_frequency,
                        disk_dissipation_age=self.disk_dissipation_age,
                        secondary_formation_age=self.disk_dissipation_age)

        secondary.configure(age=self.disk_dissipation_age,
                            companion_mass=primary.mass,
                            semimajor=binary.semimajor(self.initial_orbital_period),
                            eccentricity=0.0,
                            locked_surface=False,
                            zero_outer_inclination=True,
                            zero_outer_periapsis=True,
                            **secondary_config)




        binary.configure(age=primary.core_formation_age(),
                         semimajor=float('nan'),
                         eccentricity=float('nan'),
                         spin_angmom=numpy.array([0.0]),
                         inclination=None,
                         periapsis=None,
                         evolution_mode='LOCKED_SURFACE_SPIN')

        return binary



    def intial_second_mom(self):

        star = self.create_star(self.secondary_mass)
        planet = self.create_planet(1.0)

        binary = self.create_binary_system(star,
                                      planet,
                                      10.0)
        binary.evolve(self.disk_dissipation_age, 1e-3, 1e-6, None)

        disk_state = binary.final_state()


        planet.delete()
        star.delete()
        binary.delete()

        return numpy.array([disk_state.envelope_angmom, disk_state.core_angmom])

    def __init__(self,
                 interpolator,
                 age,
                 primary_mass,
                 secondary_mass,
                 disk_lock_frequency,
                 initial_orbital_period,
                 feh,
                 logQ):


        self.interpolator = interpolator

        self.age = age
        self.primary_mass = primary_mass
        self.secondary_mass = secondary_mass
        self.disk_lock_frequency = disk_lock_frequency
        self.initial_orbital_period = initial_orbital_period
        self.feh = feh
        self.convective_phase_lag = phase_lag(logQ)

        self.inclination = scipy.pi/2
        self.disk_dissipation_age = 5e-3
        self.wind = True
        self.planet_formation_age = 5e-3
        self.wind_saturation_frequency = 2.54
        self.diff_rot_coupling_timescale = 5e-3
        self.wind_strength = 0.17


    def __call__(self):

        print('convective_phase_lag = ', self.convective_phase_lag )

        primary = self.create_star(self.primary_mass)
        secondary = self.create_star(self.secondary_mass )



        target = Structure(age=self.age,
                           Porb=self.initial_orbital_period,  # current Porb to match
                           Wdisk=self.disk_lock_frequency,  # initial disk locking frequency
                           planet_formation_age=self.planet_formation_age)


        print('initial_second_mom = ', self.intial_second_mom())

        find_ic =InitialConditionSolver(disk_dissipation_age=self.disk_dissipation_age,
                                         evolution_max_time_step=1e-3,
                                         secondary_angmom=self.intial_second_mom(),
                                         is_secondary_star=True)

        initial_porb, final_Psurf = find_ic(target=target,
                                              primary=primary,
                                              secondary=secondary)

        primary.delete()
        secondary.delete()

        return final_Psurf



if __name__ == '__main__':

    serialized_dir = "/home/kpenev/projects/git/poet/stellar_evolution_interpolators"
    #serialized_dir = "/Users/ruskinpatel/Desktop/Research/poet/stellar_evolution_interpolators"
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        b"eccentricity_expansion_coef.txt"
    )

    age = 4.833703804484614
    primary_mass = 0.9191778701052883
    secondary_mass = 0.7315499433056569
    disk_lock_frequency = 4.554410478058684
    initial_orbital_period = 5.2663825 #this is the final orbital period quoted in lit.
    feh = -0.12783693672398547
    logQ = 8.241010816020967

    ics = ics_check(interpolator,
                    age,
                    primary_mass,
                    secondary_mass,
                    disk_lock_frequency,
                    initial_orbital_period,
                    feh,
                    logQ)

    print(ics())
