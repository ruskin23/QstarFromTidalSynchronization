import random
import scipy
from scipy.stats import norm



class struct :

    def __init__(self,**kwargs):

        for name,value in kwargs.items():
            setattr(self,name,value)



class test:

    def disp(self):

        print (self.set)

    def ini(self):

        for i in range(10):

            self.set.append(i)


    def __init__(self,observables,proposed_step,set = None):


        self.observables = observables
        self.proposed_step = proposed_step
        if set is None:

            self.set = []




    observation_data = struct(
                        Teff = struct( value = 1, sigma = 0.1 ),
                        feh = struct(value=1, sigma=0.1),
                        rvk = struct(value=1, sigma=0.1),
                        inclination = struct(value=1, sigma=0.1)
    )

    fixed_parameters = struct(
                                disk_dissipatoin_age = 5e-3
    )

    proposed_step = struct(
                        Teff_proposed_step = 0.2,
                        feh_proposed_step = 0.2,
                        rvk_proposed_step = 0.3,
                        inclination_step = 0.2
    )





class test2:

    def __init__(self,p):

        for name,value in p.items():
            setattr(self,name,value)

    def printing(self):


        print (self.b['s'])



d = dict(a = dict(v = 1, s = 0.1),
         b = dict(v = 2, s = 0.2),
         c = dict(v = 3, s = 0.3))



p = dict(psa = 0.2,
         psb = 0.3,
         psc = 0.4
         )

instance = test2(d)
#instance.printing()
n = dict()

check = p

observation_data = dict(
    Teff=dict(value=1, sigma=0.1),
    feh=dict(value=1, sigma=0.2),
    rvk=dict(value=4, sigma=0.3),
    inclination=dict(value=1, sigma=0.1)
)

#for key,value in observation_data.items():

    #print (value['value'],value['sigma'])

#for (kd,vd),(kp,vp) in zip(d.items(),p.items()):
    #print (vd['v'])
    #print (vp)
#    n[ kd] = vp*4

#for x in check:
#        print (x)#print (d[x]['v'])

width = scipy.linspace(0,1,100)

dist = norm.pdf(width,loc=0.5,scale=0.1)

class test2:

    def __init__(self,x,y):

        self.a = x
        self.b = y

    def __call__(self):

        z = self.a*self.b

        return z


arr = []

for i in range(10):
    arr.append(i)

with open('test_3.txt','w') as file:
    for x in arr:
        file.write('%s ' %x)


from mass_calculations import DeriveSecondaryMass
from inital_condition_solver import  InitialConditionSolver



____-------_______________
__________________________
__________________________



#!/usr/bin/env python3

import matplotlib

matplotlib.use('TkAgg')

# import sys
# sys.path.append('.../poet/PythonPackage')
# sys.path.append('.../poet/scripts')
import sys

sys.path.append('/Users/ruskinpatel/Desktop/Research/poet/PythonPackage')
sys.path.append('/Users/ruskinpatel/Desktop/Research/poet/scripts')

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

        print('Trying P0 = %s, Pdisk = %s'
              %
              (repr(initial_orbital_period), repr(disk_period)))

        if hasattr(self, 'binary'): self.binary.delete()

        self.secondary.select_interpolation_region(self.disk_dissipation_age)
        spin_angmom = self.secondary_angmom
        inclination = numpy.array([0.0])
        periapsis = numpy.array([0.0])
        secondary_formation_age = self.disk_dissipation_age


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

        self.binary.configure(self.primary.core_formation_age(),
                              float('nan'),
                              float('nan'),
                              numpy.array([0.0]),
                              None,
                              None,
                              'LOCKED_SURFACE_SPIN')

        self.binary.primary.detect_stellar_wind_saturation()

        self.binary.secondary.configure(
            self.target.planet_formation_age,
            self.binary.primary.mass,
            self.binary.semimajor(initial_orbital_period),
            0.0,
            spin_angmom,
            inclination,
            periapsis,
            False,
            True,
            True
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
        print('Got Porb = %s, P* = %s'
              %
              (repr(orbital_period), repr(stellar_spin_period)))
        if (numpy.isnan(orbital_period)): orbital_period = 0.0
        return orbital_period, stellar_spin_period

    def _find_porb_range(self, guess_porb_initial, disk_period):
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

        while porb_error * guess_porb_error > 0 and porb_initial < 100.0:
            if porb_error < 0:
                porb_min = porb_initial
            else:
                porb_max = porb_initial
            porb_initial *= step
            porb, psurf = self._try_initial_conditions(porb_initial,
                                                       disk_period)
            if not numpy.isnan(porb):
                porb_error = porb - self.target.Porb

        if numpy.isnan(porb_error): return numpy.nan, numpy.nan

        if porb_error < 0:
            porb_min = porb_initial
        else:
            porb_max = porb_initial
            if porb_error == 0: porb_min = porb_initial

        print('For Pdisk = %s, orbital period range: %s < Porb < %s'
              %
              (repr(disk_period), repr(porb_min), repr(porb_max)))
        return porb_min, porb_max

    def __init__(self,
                 disk_dissipation_age=None,
                 evolution_max_time_step=1.0,
                 evolution_precision=1e-6,
                 orbital_period_tolerance=1e-6,
                 spin_tolerance=1e-6,
                 secondary_angmom=None,
                 ):
        """
        Initialize the object.

        Args:
            - disk_dissipation_age:
                The age at which the disk dissipates in Gyrs.

            - evolution_max_time_step:
                The maximum timestep the evolution is allowed to make.

            - evolution_precision:
                The precision to require of the evolution.

            -secondary_angmom:
                initial angular momentum of secondary star's envelope calculated from planet-star evolution

            -

        Returns: None.
        """

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
                      return_difference=False):
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

        if numpy.isnan(porb_min):
            assert (numpy.isnan(porb_max))
            return (numpy.nan if return_difference else (numpy.nan,
                                                         numpy.nan,
                                                         numpy.nan))
        porb_initial = brentq(
            lambda porb_initial: self._try_initial_conditions(
                porb_initial,
                disk_period,
            )[0] - self.target.Porb,
            porb_min,
            porb_max,
            xtol=self.orbital_period_tolerance,
            rtol=self.orbital_period_tolerance
        )

        porb_final, spin_period = self._try_initial_conditions(
            porb_initial,
            disk_period,
        )

        spin_frequency = 2.0 * pi / spin_period

        if not return_difference:
            return spin_frequency, porb_initial, porb_final

        result = spin_frequency - 2.0 * pi / self.target.Psurf
        if (
                abs(result)
                <
                abs(self._best_initial_conditions.spin_frequency
                    -
                    2.0 * pi / self.target.Psurf)
        ):
            self._best_initial_conditions.spin_frequency = spin_frequency
            self._best_initial_conditions.orbital_period = porb_final
            self._best_initial_conditions.initial_orbital_period = porb_initial
            self._best_initial_conditions.disk_period = disk_period
        return result

    def __call__(self, target, primary, secondary):
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
                    - Wdisk:
                        the disk locking frequency.
                The following optional attributes can be specified:
                    - disk_dissipation_age:
                        The age at which the disk dissipates in Gyrs. If not
                        specified, it must have been defined when this solver
                        was initialized.
            - primary:
                should be instance of evolve_interface.EvolvingStar and its dissipative properties
                should be defined.

            - secondary:
                should be instance of evolve_interface.EvolvingStar and its dissipative properties
                should be defined.

        Returns:
            - porb_initial:
                Initial orbital period.

            - pstar:
                stellar spin.

            Further, the solver object has an attribute named 'binary' (an
            instance of (evolve_interface.Binary) which was evolved from
            the initial conditions found to most closely reproduce the
            specified target configuration.
        """


        self.target = target
        self.primary = primary
        self.secondary = secondary

        if not hasattr(self.target, 'disk_dissipation_age'):
            self.target.disk_dissipation_age = self.disk_dissipation_age

        Wdisk = (target.Wdisk if hasattr(target, 'Wdisk')
                     else 2.0 * pi / target.Pdisk)
        Wstar, Porb_initial, Porb_now = self.stellar_wsurf(Wdisk,target.Porb)
        return Porb_initial, 2.0 * pi / Wstar

