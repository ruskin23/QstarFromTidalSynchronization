#!/usr/bin/env python3

import pickle

import sys

from pathlib import Path
home_dir=str(Path.home())

if home_dir=='/home/rxp163130':
    poet_path=home_dir+'/poet'
if home_dir=='/home/ruskin':
    poet_path=home_dir+'/projects/poet'
if home_dir=='/home1/06850/rpatel23':
    work_dir='/work/06850/rpatel23/stampede2'
    poet_path=work_dir+'/poet'

sys.path.append(poet_path+'/PythonPackage')
sys.path.append(poet_path+'/scripts')


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
import scipy

class InitialConditionSolver:
    """Find initial conditions which reproduce a given system now."""

    def _try_initial_conditions(self, initial_orbital_period, disk_period):
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

        print('\nTrying P0 = %s, Pdisk = %s' %(repr(initial_orbital_period), repr(disk_period)))

        if hasattr(self, 'binary'): self.binary.delete()

        if self.is_secondary_star is True:
            self.secondary.select_interpolation_region(self.disk_dissipation_age)
            spin_angmom = self.secondary_angmom
            inclination = numpy.array([0.0])
            periapsis = numpy.array([0.0])
            secondary_formation_age = self.disk_dissipation_age
            secondary_config_age = self.disk_dissipation_age


        else:
            spin_angmom = numpy.array([0.0])
            inclination = None
            periapsis = None
            secondary_formation_age = self.target.planet_formation_age

        self.binary = Binary(
            primary=self.primary,
            secondary=self.secondary,
            initial_orbital_period=initial_orbital_period,
            initial_eccentricity=0.0,
            initial_inclination=0.0,
            disk_lock_frequency=2.0 * numpy.pi / disk_period,
            disk_dissipation_age=self.disk_dissipation_age,
            secondary_formation_age=secondary_formation_age
        )

        self.binary.primary.select_interpolation_region(self.primary.core_formation_age())
        if self.is_secondary_star is True: self.binary.secondary.detect_stellar_wind_saturation()


        self.binary.configure(  age =  self.primary.core_formation_age(),
                                semimajor=float('nan'),
                                eccentricity=float('nan'),
                                spin_angmom=numpy.array([0.0]),
                                inclination=None,
                                periapsis=None,
                                evolution_mode='LOCKED_SURFACE_SPIN'

                                )
        self.binary.primary.detect_stellar_wind_saturation()

        self.binary.secondary.configure(
            age=self.disk_dissipation_age,
            companion_mass=self.binary.primary.mass,
            semimajor=self.binary.semimajor(initial_orbital_period),
            eccentricity=0.0,
            spin_angmom=spin_angmom,
            inclination=inclination,
            periapsis=periapsis,
            locked_surface=False,
            zero_outer_inclination=True,
            zero_outer_periapsis=True
        )

        print ("BINARY CONFIGURATION COMPLETE")

        self.binary.evolve(
            self.target.age,
            self.evolution_max_time_step,
            self.evolution_precision,
            None
        )

        print ("BINARY EVOLUTION COMPLETE")

        final_state = self.binary.final_state()

        assert (final_state.age == self.target.age)

        orbital_period = self.binary.orbital_period(final_state.semimajor)

        stellar_spin_period = (
                2.0 * pi
                *
                self.binary.primary.envelope_inertia(final_state.age)
                /
                final_state.primary_envelope_angmom
        )
        print('Got Porb = %s, P* = %s'     % (repr(orbital_period), repr(stellar_spin_period)))
        if (numpy.isnan(orbital_period)): orbital_period = 0.0
        return orbital_period, stellar_spin_period

    def _find_porb_range(self, guess_porb_initial, disk_period):

        porb_min, porb_max = numpy.nan, numpy.nan
        porb_initial = guess_porb_initial
        try:
            porb, psurf = self._try_initial_conditions(porb_initial, disk_period)
        except AssertionError:
            return scipy.nan,scipy.nan
        porb_error = porb - self.target.Porb
        guess_porb_error = porb_error
        step = 2.0 if guess_porb_error < 0 else 0.5

        while porb_error * guess_porb_error > 0 and porb_initial < 100.0:
            if porb_error < 0:
                porb_min = porb_initial
            else:
                porb_max = porb_initial
            porb_initial *= step
            try:
                porb, psurf = self._try_initial_conditions(porb_initial,
                                                       disk_period)
                if not numpy.isnan(porb):
                    porb_error = porb - self.target.Porb
            except AssertionError:
                return scipy.nan,scipy.nan

        if numpy.isnan(porb_error): return numpy.nan, numpy.nan

        if porb_error < 0:
            porb_min = porb_initial
        else:
            porb_max = porb_initial
            if porb_error == 0: porb_min = porb_initial

        print('\nFOR Pdisk = %s, ORBITAL PERIOD RANGE: %s < Porb < %s'       % (repr(disk_period), repr(porb_min), repr(porb_max)))
        return porb_min, porb_max

    def __init__(self,
                 planet_formation_age=None,
                 disk_dissipation_age=None,
                 evolution_max_time_step=None,
                 evolution_precision=1e-6,
                 secondary_angmom=None,
                 is_secondary_star=None):

        self.disk_dissipation_age = disk_dissipation_age
        self.evolution_max_time_step = evolution_max_time_step
        self.evolution_precision = evolution_precision
        self.secondary_angmom = secondary_angmom
        self.is_secondary_star = is_secondary_star


    def stellar_period(self,
                      wdisk,
                      orbital_period_guess):

        disk_period = 2.0 * pi / wdisk

        porb_min, porb_max = self._find_porb_range(orbital_period_guess,
                                                   disk_period)


        if numpy.isnan(porb_min):
            assert (numpy.isnan(porb_max))
            return numpy.nan,numpy.nan

        try:
            porb_initial = brentq(
                lambda porb_initial: self._try_initial_conditions(
                    porb_initial,
                    disk_period,
                )[0] - self.target.Porb,
                porb_min,
                porb_max
            )
        except ValueError:
            return scipy.nan.scipy.nan

        porb_final, spin_period = self._try_initial_conditions(
            porb_initial,
            disk_period,
        )

        return spin_period,porb_initial

    def __call__(self, target, primary, secondary):

        self.target = target
        self.primary = primary
        self.secondary = secondary

        Wdisk = target.Wdisk
        spin_period, Porb_initial= self.stellar_period(Wdisk,target.Porb)

        print ("\n Results Spin_period = %s, Porb_initial = %s, Porb_now = %s" % (repr(spin_period),repr(Porb_initial),repr(target.Porb)))

        return spin_period


