#!/usr/bin/env python3

import pickle

import sys

sys.path.append('/home/kpenev/projects/git/poet/PythonPackage')
sys.path.append('/home/kpenev/projects/git/poet/scripts')


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
        #if hasattr(self, 'binary'): self.binary.delete()
        if initial_condition[1]>0.45 or initial_condition[1]<0:
            print('Cannot accept eccentricity>0.45')
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
        assert (self.final_state.age == self.target.age)

        self.orbital_period = self.binary.orbital_period(self.final_state.semimajor)
        self.eccentricity = self.final_state.eccentricity


        if (numpy.isnan(self.orbital_period)): self.orbital_period = 0.0

        self.delta_p = self.orbital_period-self.target.Porb
        self.delta_e = self.eccentricity-self.target.eccentricity

        print(self.delta_p,self.delta_e)

        self.spin =  (
                2.0 * pi
                *
                self.binary.primary.envelope_inertia(self.final_state.age)
                /
                self.final_state.primary_envelope_angmom
        )

        self.binary.delete()
        print(self.spin)
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

        ecc=numpy.linspace(e,0.4,5)
        porb=numpy.linspace(p,p*2.0,5)

        print('eccentricity = ', ecc)
        print('Porb = ', porb)

        with open('porb_e.txt','w',1) as f:
            for ex in ecc:
                for px in porb:
                    print(px)
                    a=[px,ex]
                    try:
                        result=self._try_initial_conditions(a)
                        f.write(repr(px)+'\t'+
                            repr(ex)+'\t'+
                            repr(self.orbital_period)+'\t'+
                            repr(self.eccentricity)+'\n')
                    except:continue
        """
        while True:
            try:
                print('solving for p and e')
                sol = optimize.root(self._try_initial_conditions,
                                    [p,e],
                                    method='lm',
                                    tol=1e-4,
                                    options={'xtol':1e-5}
                                    )
                sol_p,sol_e=sol.x
                break

            except Exception:
                if self.delta_p<1e-5 and self.delta_e<1e-5:
                    sol_p,sol_e=self.p_initial,self.e_inital
                    break
                else:
                    if self.e_inital:e=self.e_inital+0.001
                    else:e=e+0.001
                    continue

        solutions=dict()

        solutions['spin']=self.spin
        solutions['Porb_inital']=sol_p
        solutions['e_initial']=sol_e
        solutions['Porb_current']=self.orbital_period
        solutions['e_current']=self.eccentricity
        solutions['delta_p']=self.delta_p
        solutions['delta_e']=self.delta_p
        """
        #return solutions
