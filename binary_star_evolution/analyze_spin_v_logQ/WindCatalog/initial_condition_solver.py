#!/usr/bin/env python3

import pickle

import sys
import os
import os.path

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
import scipy
from scipy import optimize
from scipy.optimize import brentq
import numpy

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
        sys.stdout.flush()
        #if hasattr(self, 'binary'): self.binary.delete()
        if initial_condition[1]>0.45 or initial_condition[1]<0:
            print('unacceptable values in initial_condition')
            sys.stdout.flush()
            return scipy.nan,scipy.nan

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
        sys.stdout.flush()

        self.binary.evolve(
                self.target.age,
                self.evolution_max_time_step,
                self.evolution_precision,
                None
            )

        self.p_initial = initial_condition[0]
        self.e_inital = initial_condition[1]

        print ("BINARY EVOLUTION COMPLETE")

        self.final_state = self.binary.final_state()
        print('Current Age = ', self.final_state.age)
        sys.stdout.flush()
        assert (self.final_state.age == self.target.age)

        self.orbital_period = self.binary.orbital_period(self.final_state.semimajor)
        self.eccentricity = self.final_state.eccentricity


        if (numpy.isnan(self.orbital_period)): self.orbital_period = 0.0

        self.delta_p = self.orbital_period-self.target.Porb
        self.delta_e = self.eccentricity-self.target.eccentricity

        print('Final orbital period = ',self.orbital_period)
        print('Final Eccentricity = ', self.eccentricity)

        print('Difference between final and target p=%s , e=%s '
              %(repr(self.delta_p),repr(self.delta_e)))


        sys.stdout.flush()
        self.spin =  (
                2.0 * pi
                *
                self.binary.primary.envelope_inertia(self.final_state.age)
                /
                self.final_state.primary_envelope_angmom
        )

        self.binary.delete()
        print(self.spin)
        sys.stdout.flush()
        return self.delta_p,self.delta_e

    def __init__(self,
                 system=None,
                 print_cfile=None,
                 breaks=None,
                 disk_dissipation_age=None,
                 evolution_max_time_step=None,
                 evolution_precision=1e-6,
                 secondary_angmom=None,
                 is_secondary_star=None):
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
        self.system=system
        self.print_cfile=print_cfile
        self.breaks=breaks
        self.disk_dissipation_age = disk_dissipation_age
        self.evolution_max_time_step = evolution_max_time_step
        self.evolution_precision = evolution_precision
        self.secondary_angmom = secondary_angmom
        self.is_secondary_star = is_secondary_star
        self.delta_p=1.0
        self.delta_e=1.0

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
                    - Pdisk | Wdisk:
                        The stellar surface spin period to reproduce or the
                        disk locking period or the disk locking frequency.
            - primary

            - secondary

        Returns:
            - porb_initial:
                Initial orbital period.

            - psurf:
                Current priamry star spin matching initial orbital period and
                initial_disk_fequency

            Further, the solver object has an attribute named 'binary' (an
            instance of (evolve_interface.Binary) which was evolved from
            the initial conditions found to most closely reproduce the
            specified target configuration.
        """


        self.target = target
        self.primary = primary
        self.secondary = secondary

        self.disk_lock_frequency = (target.Wdisk if hasattr(target, 'Wdisk')
                                    else 2*pi/target.Pdisk)


        e=target.eccentricity
        p=target.Porb

        solutions=dict()

        while True:
            try:
                print('solving for p and e')
                sol = optimize.root(self._try_initial_conditions,
                                    [p,e],
                                    method='lm',
                                    options={'xtol':1e-5,
                                             'ftol':1e-5}
                                    )
                sol_p,sol_e=sol.x
                break

            except Exception as err:
                print('err = ', err)
                sys.stdout.flush()
                #break
                return solutions

        solutions['spin']=self.spin
        solutions['Porb_inital']=sol_p
        solutions['e_initial']=sol_e
        solutions['Porb_current']=self.orbital_period
        solutions['e_current']=self.eccentricity
        solutions['delta_p']=self.delta_p
        solutions['delta_e']=self.delta_e

        return solutions
