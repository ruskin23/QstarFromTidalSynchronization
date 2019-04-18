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
from scipy import optimize
from scipy.optimize import brentq
import numpy
import matplotlib.pyplot as plt

wsun = 0.24795522138

class InitialConditionSolver:
    """Find initial conditions which reproduce a given system now."""

    def try_initial_conditions(self, initial_condition):
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

        print ("BINARY EVOLUTION COMPLETE")

        self.final_state = self.binary.final_state()
        assert (self.final_state.age == self.target.age)

        self.evolution=self.binary.get_evolution()
        self.evolution_age=self.evolution.age
        self.wenv_primary=(self.evolution.primary_envelope_angmom / self.binary.primary.envelope_inertia(self.evolution.age)) / wsun

        print('final semimajor = ', self.final_state.semimajor)
        print('final age = ', self.final_state.age)
        print('final primary envelope angular momentum = ', self.final_state.primary_envelope_angmom)
        print('final primary core angular momentum', self.final_state.primary_core_angmom)
        print('final secondary envelope angular momentum = ', self.final_state.secondary_envelope_angmom)
        print('final secondary core angular momentum', self.final_state.secondary_core_angmom)


        self.orbital_period = self.binary.orbital_period(self.final_state.semimajor)
        self.eccentricity = self.final_state.eccentricity

        if (numpy.isnan(self.orbital_period)): self.orbital_period = 0.0

        self.delta_p = self.orbital_period-self.target.Porb
        self.detla_e = self.eccentricity-self.target.eccentricity

        print(self.delta_p,self.detla_e)

        self.spin =  (
                2.0 * pi
                *
                self.binary.primary.envelope_inertia(self.final_state.age)
                /
                self.final_state.primary_envelope_angmom
        )

        return self.delta_p,self.detla_e

    def __init__(self,
                 planet_formation_age=None,
                 disk_dissipation_age=None,
                 evolution_max_time_step=None,
                 evolution_precision=1e-6,
                 secondary_angmom=None,
                 is_secondary_star=None):

        self.final_state=None
        self.disk_dissipation_age = disk_dissipation_age
        self.evolution_max_time_step = evolution_max_time_step
        self.evolution_precision = evolution_precision
        self.secondary_angmom = secondary_angmom
        self.is_secondary_star = is_secondary_star

        self.evolution = None
        self.wenv_primary = None
        self.evolution_age = None

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


        check_intial_condition=[target.Porb,target.eccentricity]

        self.try_initial_conditions(check_intial_condition)
        #plt.semilogx(self.evolution_age,self.wenv_primary,'-r')
        age_1 = self.evolution_age
        primary_envelope_angmom_1 = numpy.array(self.wenv_primary)

        self.try_initial_conditions(check_intial_condition)
        #plt.semilogx(self.evolution_age,self.wenv_primary,'-g')
        age_2 = self.evolution_age
        primary_envelope_angmom_2 = numpy.array(self.wenv_primary)

        print(len(primary_envelope_angmom_1))
        print(len(primary_envelope_angmom_2))

        print(len(age_1))
        print(len(age_2))

        with open('evolve_1.txt','w') as f:
            i=1
            for age,env in zip(age_1,primary_envelope_angmom_1):
                f.write(repr(i) + '\t' + repr(age) + '\t' + repr(env) + '\n')
                i=i+1

        with open('evolve_2.txt','w') as f:
            i=1
            for age,env in zip(age_2,primary_envelope_angmom_2):
                f.write(repr(i) + '\t' + repr(age) + '\t' + repr(env) + '\n')
                i=i+1

        #difference_array = primary_envelope_angmom_1-primary_envelope_angmom_2


        #plt.plot(self.evolution_age,difference_array)
        #plt.show()

