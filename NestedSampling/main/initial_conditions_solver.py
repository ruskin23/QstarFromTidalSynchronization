#!/usr/bin/env python3 -u

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
from scipy import optimize
from scipy.optimize import brentq
import numpy
import scipy

class InitialConditionSolver:
    """Find initial conditions which reproduce a given system now."""

    def _try_initial_conditions(self, initial_condition):
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

        #print('\nTrying P0 = %s, Pdisk = %s' %(repr(initial_orbital_period), repr(disk_period)))
        print('\nTrying Porb_initial = %s, e_initial =%s'
              %(repr(initial_condition[0]), repr(initial_condition[1])))
        if initial_condition[1]>0.45 or initial_condition[1]<0 or initial_condition[0]<0:
            print('Cannot accept eccentricity > 0.45')
            return scipy.nan, scipy.nan

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
            initial_orbital_period=initial_condition[0],
            initial_eccentricity=initial_condition[1],
            initial_inclination=0.0,
            disk_lock_frequency=self.disk_lock_frequency,
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
            semimajor=self.binary.semimajor(initial_condition[0]),
            eccentricity=initial_condition[1],
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
            None,
            timeout=600
        )

        print ("BINARY EVOLUTION COMPLETE")

        self.p_initial=initial_condition[0]
        self.e_initial=initial_condition[1]

        self.final_state = self.binary.final_state()
        assert (self.final_state.age == self.target.age)

        self.orbital_period = self.binary.orbital_period(self.final_state.semimajor)
        self.eccentricity = self.final_state.eccentricity



        self.delta_p = self.orbital_period-self.target.Porb
        self.delta_e = self.eccentricity-self.target.eccentricity

        if numpy.logical_or(numpy.isnan(self.delta_p),(numpy.isnan(self.delta_e)))==True:
            print('Binary system was destroyed')
            raise ValueError

        print(self.delta_p,self.delta_e)

        self.spin =  (
                2.0 * pi
                *
                self.binary.primary.envelope_inertia(self.final_state.age)
                /
                self.final_state.primary_envelope_angmom
        )

        print(self.spin)

        self.binary.delete()
        sys.stdout.flush()

        return self.delta_p,self.delta_e

    def __init__(self,
                 planet_formation_age=None,
                 disk_dissipation_age=None,
                 evolution_max_time_step=None,
                 evolution_precision=1e-6,
                 secondary_angmom=None,
                 is_secondary_star=None):
        """
        Initialize the object.
        """

        self.disk_dissipation_age = disk_dissipation_age
        self.evolution_max_time_step = evolution_max_time_step
        self.evolution_precision = evolution_precision
        self.secondary_angmom = secondary_angmom
        self.is_secondary_star = is_secondary_star
        self.e_initial = 0
        self.p_initial = 0
        self.initial_orbital_period_sol=None
        self.initial_eccentricity_sol=None
        self.bad_solution=False
        self.gsl_flag=False

    def __call__(self, target, primary, secondary):
        """
        Find initial conditions which reproduce the given system now.

        """


        self.target = target
        self.primary = primary
        self.secondary = secondary

        self.disk_lock_frequency = (target.Wdisk if hasattr(target, 'Wdisk')
                                    else 2*pi/target.Pdisk)


        e=self.target.eccentricity
        p=self.target.Porb
        while True:
            try:
                sol = optimize.root(
                                self._try_initial_conditions,
                                [p,e],
                                method='lm'
                                )
                break
            except:
                self.spin=scipy.nan
                self.delta_e=scipy.nan
                self.delta_p=scipy.nan
                self.orbital_period=p
                self.eccentricity=e
                break


        print(self.spin)

        if numpy.isnan(self.spin):
            self.initial_orbital_period_sol=scipy.nan
            self.initial_eccentricity_sol=scipy.nan

        else:
            self.initial_orbital_period_sol,self.initial_eccentricity_sol = sol.x

        if abs(self.delta_p)>0.1 or abs(self.delta_e)>0.1:
            self.bad_solution=True
            self.spin=scipy.nan

        print('spin from ICS = ', self.spin)
        return self.spin





