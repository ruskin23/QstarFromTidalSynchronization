#!/usr/bin/env python3

import matplotlib

matplotlib.use('TkAgg')

import sys
sys.path.append('.../poet/PythonPackage')
sys.path.append('.../poet/scripts')
#import sys
#sys.path.append('/Users/ruskinpatel/Desktop/Research/poet/PythonPackage')
#sys.path.append('/Users/ruskinpatel/Desktop/Research/poet/scripts')



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

        if self.is_secondary_star is True:
            self.secondary.select_interpolation_region(self.disk_dissipation_age)
            spin_angmom=self.secondary_angmom
            inclination=numpy.array([0.0])
            periapsis=numpy.array([0.0])
            secondary_formation_age = self.disk_dissipation_age
            secondary_config_age = self.disk_dissipation_age


        else:
            spin_angmom=numpy.array([0.0])
            inclination=None
            periapsis=None
            secondary_formation_age = self.target.planet_formation_age





        self.binary = Binary(
            primary = self.primary,
            secondary = self.secondary,
            initial_orbital_period = initial_orbital_period,
            initial_eccentricity = 0.0,
            initial_inclination = 0.0,
            disk_lock_frequency = 2.0 * numpy.pi / disk_period,
            disk_dissipation_age = self.disk_dissipation_age,
            secondary_formation_age =secondary_formation_age
        )

        self.binary.primary.select_interpolation_region(self.primary.core_formation_age())


        self.binary.configure(  age =  self.primary.core_formation_age(),
                                semimajor=float('nan'),
                                eccentricity=float('nan'),
                                spin_angmom=numpy.array([0.0]),
                                inclination = None,
                                periapsis = None,
                                evolution_mode='LOCKED_SURFACE_SPIN'

                                )
                            

        self.binary.primary.detect_stellar_wind_saturation()
        if self.is_secondary_star is  True:  self.binary.secondary.detect_stellar_wind_saturation()


        self.binary.secondary.configure(
            age = self.target.planet_formation_age,
            companion_mass = self.binary.primary.mass,
            semimajor = self.binary.semimajor(initial_orbital_period),
            eccentricity = 0.0,
            spin_angmom = spin_angmom,
            inclination = inclination,
            periapsis = periapsis,
            locked_surface = False,
            zero_outer_inclination = True,
            zero_outer_periapsis = True
        
                 ) 

        print ("BINARY CONFIGURATION COMPLETE")




        self.binary.evolve(
            self.target.age,
            self.evolution_max_time_step,
            self.evolution_precision,
            None,
            True
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
            final_state.primary_envelope_angmom
        )
        print('Got Porb = %s, P* = %s'
              %
              (repr(orbital_period), repr(stellar_spin_period)))
        if(numpy.isnan(orbital_period)) : orbital_period = 0.0
        return orbital_period, stellar_spin_period

    def _find_porb_range(self, guess_porb_initial, disk_period) :
        """
        Find initial orbital period range where final porb error flips sign.

        Args:
            - guess_porb_initial:
                An initial guess for where the sign change occurs.

            - disk_period:
                The disk locking period to assume during the search.

        Returns:
            A pair of initial orbital periods for which the sign of the final
            orbital period error changes.
        """

        porb_min, porb_max = numpy.nan, numpy.nan
        porb_initial = guess_porb_initial
        porb, psurf = self._try_initial_conditions(porb_initial, disk_period)
        porb_error = porb - self.target.Porb
        guess_porb_error = porb_error
        step = 2.0 if guess_porb_error < 0 else 0.5

        while porb_error * guess_porb_error > 0 and porb_initial < 100.0 :
            if porb_error < 0 : porb_min = porb_initial
            else : porb_max = porb_initial
            porb_initial *= step
            porb, psurf = self._try_initial_conditions(porb_initial,
                                                       disk_period)
            if not numpy.isnan(porb) :
                porb_error = porb - self.target.Porb

        if numpy.isnan(porb_error) : return numpy.nan, numpy.nan

        if porb_error < 0 : porb_min = porb_initial
        else :
            porb_max = porb_initial
            if porb_error == 0 : porb_min = porb_initial

        print('For Pdisk = %s, orbital period range: %s < Porb < %s'
              %
              (repr(disk_period), repr(porb_min), repr(porb_max)))
        return porb_min, porb_max

    def __init__(self,
                 planet_formation_age = None,
                 disk_dissipation_age = None,
                 evolution_max_time_step = 1.0,
                 evolution_precision = 1e-4,
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

    def stellar_wsurf(self,
                      wdisk,
                      orbital_period_guess,
                      return_difference=False) :
        """
        The stellar spin frquency when reproducing current porb.

        Args:
            - disk_frequency: The angular velocity of the star when it forms.
            - orbital_period_guess: A best guess value for the initial
                                    orbital period.
            - return_difference: If True, instead of the actual stellar
                                 angulare velocity, the function returns the
                                 difference from the observed value.

        Returns:
            - spin_frequency:
                The angular velocity with which the star spins at the present
                age for an evolution scenario which reproduces the current
                orbital period. Or the difference between the spin frequency
                and the target spin frequency if return_difference is True.

            The following are returned only if return_difference is False:
            - porb_initial:
                The initial orbital period which reproduces the specified
                final orbital period as close as possible.

            - porb_final:
                The closest final orbital period found (starting with
                porb_initial).
        """

        disk_period = 2.0 * pi / wdisk

        porb_min, porb_max = self._find_porb_range(orbital_period_guess,
                                                   disk_period)

        if numpy.isnan(porb_min) :
            assert(numpy.isnan(porb_max))
            return (numpy.nan if return_difference else (numpy.nan,
                                                         numpy.nan,
                                                         numpy.nan))
        porb_initial = brentq(
            lambda porb_initial : self._try_initial_conditions(
                porb_initial,
                disk_period,
            )[0] - self.target.Porb,
            porb_min,
            porb_max,
            xtol = self.orbital_period_tolerance,
            rtol = self.orbital_period_tolerance
        )

        porb_final, spin_period = self._try_initial_conditions(
            porb_initial,
            disk_period,
        )

        spin_frequency = 2.0 * pi / spin_period

        if not return_difference :
            return spin_frequency, porb_initial, porb_final

        result = spin_frequency - 2.0 * pi / self.target.Psurf
        if (
                abs(result)
                <
                abs(self._best_initial_conditions.spin_frequency
                    -
                    2.0 * pi / self.target.Psurf)
        ) :
            self._best_initial_conditions.spin_frequency = spin_frequency
            self._best_initial_conditions.orbital_period = porb_final
            self._best_initial_conditions.initial_orbital_period = porb_initial
            self._best_initial_conditions.disk_period = disk_period
        return result

    def __call__(self, target,primary , secondary):
        """
        Find initial conditions which reproduce the given system now.

        Args:
            - target:
                The target configuration to reproduce by tuning the the
                initial conditions for.
                The following attributes must be defined:
                    - age:
                        The age at which the system configuration is known.
                    - Porb:
                        The orbital period to reproduce.
                    - Psurf | Pdisk | Wdisk:
                        The stellar surface spin period to reproduce or the
                        disk locking period or the disk locking frequency.
                The following optional attributes can be specified:
                    - planet_formation_age:
                        The age at which the planet forms in Gyrs. If not
                        specified the planet is assumed to form either
                        '-past_lifetime' Gyrs before '-age' or as soon as
                        the disk dissipates.
                    - past_lifetime:
                        An alternative way of specifying when the planet
                        forms. If the '-planet_formation-age' attribute is
                        not defined, the planet is assumed to form this many
                        Gyr before '-age'.
                    - disk_dissipation_age:
                        The age at which the disk dissipates in Gyrs. If not
                        specified, it must have been defined when this solver
                        was initialized.
            - star:
                The star to use in the evolution, should be instance of
                evolve_interface.EvolvingStar and its dissipative properties
                should be defined.

            - planet:
                The planet to use in the evolution. Should be an instance of
                evolve_interface.LockedPlanet

        Returns:
            - porb_initial:
                Initial orbital period.

            - pdisk | pstar:
                Initial disk period if matching observed stellar spin or
                stellar spin if initial disk period is specified.

            Further, the solver object has an attribute named 'binary' (an
            instance of (evolve_interface.Binary) which was evolved from
            the initial conditions found to most closely reproduce the
            specified target configuration.
        """

        def reset_best_initial_conditions():
            """Reset the entries in self._best_initial_conditions."""

            self._best_initial_conditions = Structure(
                spin_frequency=numpy.inf,
                orbital_period=numpy.nan,
                initial_orbital_period=numpy.nan,
                disk_period=numpy.nan
            )

        def get_initial_grid():
            """Tabulate stellar spin errors for a grid of disk periods."""

            reset_best_initial_conditions()

            wdisk_grid = [-200.0 * pi,
                          -20.0 * pi,
                          -2.0 * pi,
                          -0.2 * pi,
                          -0.02 * pi,
                          -0.002 * pi,
                          0.002 * pi,
                          0.02 * pi,
                          0.2 * pi,
                          2.0 * pi,
                          20.0 * pi,
                          200.0 * pi]

            stellar_wsurf_residual_grid = [
                self.stellar_wsurf(wdisk, target.Porb, True)
                for wdisk in wdisk_grid
            ]

            print('##    %25s %25s\n' % ('disk_period',
                                         'current_stellar_spin'))
            for wdisk, wsurf_residual in zip(wdisk_grid,
                                             stellar_wsurf_residual_grid):
                print('##    %25.16e %25.16e\n'
                      %
                      (2.0 * pi / wdisk,
                       wsurf_residual + 2.0 * pi / self.target.Psurf))
            print('## Target current stellar spin: '
                  +
                  repr(2.0 * pi / self.target.Psurf)
                  +
                  '\n')

            return wdisk_grid, stellar_wsurf_residual_grid

        self.target = target
        self.primary = primary
        self.secondary = secondary

        if not hasattr(self.target, 'disk_dissipation_age'):
            self.target.disk_dissipation_age = self.disk_dissipation_age
        if not hasattr(self.target, 'planet_formation_age'):
            if hasattr(self.target, 'past_lifetime'):
                self.target.planet_formation_age = (self.target.age
                                                    -
                                                    self.target.past_lifetime)
            else:
                self.target.planet_formation_age = self.planet_formation_age

        if not hasattr(target, 'Psurf'):
            Wdisk = (target.Wdisk if hasattr(target, 'Wdisk')
                     else 2.0 * pi / target.Pdisk)
            Wstar, Porb_initial, Porb_now = self.stellar_wsurf(Wdisk,
                                                               target.Porb)
            return Porb_initial, 2.0 * pi / Wstar
        else:
            wdisk_grid, stellar_wsurf_residual_grid = get_initial_grid()

            nsolutions = 0
            for i in range(len(wdisk_grid) - 1):
                if (
                        stellar_wsurf_residual_grid[i]
                        *
                        stellar_wsurf_residual_grid[i + 1]
                        <
                        0
                ):
                    wdisk = brentq(
                        f=self.stellar_wsurf,
                        a=wdisk_grid[i],
                        b=wdisk_grid[i + 1],
                        args=(target.Porb,
                              True),
                        xtol=self.spin_tolerance,
                        rtol=self.spin_tolerance
                    )

                    nsolutions += 1

            return (self._best_initial_conditions.initial_orbital_period,
                    self._best_initial_conditions.disk_period)





def create_planet(mass=(constants.M_jup / constants.M_sun).to('')):
    """Return a configured planet to use in the evolution."""

    planet = LockedPlanet( mass=mass,radius=(constants.R_jup / constants.R_sun).to('') )
    return planet


def create_star(mass, interpolator, convective_phase_lag, wind=True):
    star = EvolvingStar(mass=mass,
                            metallicity=  -0.06,
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

    return star


def create_binary_system(primary,
                         secondary,
                         disk_lock_frequency,
                         initial_orbital_period,
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


    if isinstance(secondary, EvolvingStar):
        secondary.detect_stellar_wind_saturation()

    primary.select_interpolation_region(primary.core_formation_age())
    primary.detect_stellar_wind_saturation()

    binary = Binary(primary=primary,
                    secondary=secondary,
                    initial_orbital_period=initial_orbital_period,
                    initial_eccentricity=0.0,
                    initial_inclination=0.0,
                    disk_lock_frequency=disk_lock_frequency,
                    disk_dissipation_age=disk_dissipation_age,
                    secondary_formation_age=disk_dissipation_age)

    secondary.configure(age=disk_dissipation_age,
                        companion_mass=primary.mass,
                        semimajor=binary.semimajor(initial_orbital_period),
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


def test_ic_solver(interpolator,convective_phase_lag,wind):
    """Find initial condition to reproduce some current state and plot."""



    tdisk = 5e-3
    age = 4.60 
    primary_mass = 1.0943124828794022
    secondary_mass = 0.8266661262872534
    initial_disk_period =1.4064967495370835
    initial_orbital_period = 5.2663825 #this is the final orbital period quoted in lit.



    star = create_star(secondary_mass, interpolator=interpolator, convective_phase_lag=convective_phase_lag, wind=wind)
    planet = create_planet(1.0)

    binary = create_binary_system(star,
                                  planet,
                                  2.0 * numpy.pi / initial_disk_period,
                                  initial_orbital_period,
                                  tdisk)

    binary.evolve(tdisk, 1e-3, 1e-6, None)

    disk_state = binary.final_state()

    print("FINISHED PLANET-STAR EVOLUTION")

    planet.delete()
    star.delete()
    binary.delete()



    primary = create_star(primary_mass, interpolator, convective_phase_lag, wind = wind)
    secondary = create_star(secondary_mass, interpolator, convective_phase_lag, wind = wind)


    find_ic = InitialConditionSolver(disk_dissipation_age=tdisk,
                                         evolution_max_time_step=1e-3,
                                         secondary_angmom=numpy.array(
                                             [disk_state.envelope_angmom, disk_state.core_angmom]),
                                         is_secondary_star=True)




    target = Structure(age=age,
                           Porb=initial_orbital_period,
                           Wdisk=2 * numpy.pi / initial_disk_period,
                           planet_formation_age=5e-3)

    initial_porb, initial_psurf = find_ic(target=target,
                                              primary=primary,
                                              secondary=secondary)
    

    primary.delete()
    secondary.delete()

    return initial_psurf

    


if __name__ == '__main__':
    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        b"eccentricity_expansion_coef.txt"
    )
    serialized_dir = '/home/kpenev/projects/git/poet/stellar_evolution_interpolators'
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    # test_evolution(interpolator, phase_lag(6.0))
    #logQ = [6.0]

    logQ = float(sys.argv[1])
    name = 'Psurf_values_' + str(logQ)
    with open(name,'w') as f:

        lag = phase_lag(logQ)
        print( "calculatingforthevalues = " + repr(lag) + "\t"+ repr(logQ) + "\n")
        Pspin = test_ic_solver(interpolator,lag,True)
        f.write(repr(lag) + "\t" + repr(logQ) + "\t" + repr(Pspin) + "\n" )



