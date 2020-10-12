#!/usr/bin/env python3 -u

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
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from orbital_evolution.binary import Binary
from orbital_evolution.transformations import phase_lag
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet
from stellar_evolution.derived_stellar_quantities import\
    TeffK,\
    LogGCGS,\
    RhoCGS
from basic_utils import Structure

import pickle
import numpy
import scipy
from astropy import units, constants


class SpinPeriod():

    def create_planet(self,mass=(constants.M_jup / constants.M_sun).to('')): 
        """Return a configured planet to use in the evolution."""

        planet = LockedPlanet(mass=mass, radius=(constants.R_jup / constants.R_sun).to(''))
        return planet

    def create_star(self, mass, dissipation):
        star = EvolvingStar(mass=mass,
                            metallicity=self.feh,
                            wind_strength=self.wind_strength if self.wind else 0.0, 
                            wind_saturation_frequency=self.wind_saturation_frequency, 
                            diff_rot_coupling_timescale=self.diff_rot_coupling_timescale, 
                            interpolator=self.interpolator)

        if dissipation == 1:
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
                             initial_orbital_period,
                             initial_eccentricity,
                             secondary_angmom=None):
        """Create a binary system to evolve from the given objects."""

        if isinstance(secondary, LockedPlanet):
            spin_angmom=numpy.array([0.0])
            inclination=None
            periapsis=None
        else:
            secondary.select_interpolation_region(self.disk_dissipation_age)
            spin_angmom=secondary_angmom
            inclination=numpy.array([0.0])
            periapsis=numpy.array([0.0])

        binary = Binary(primary=primary,
                        secondary=secondary,
                        initial_orbital_period=initial_orbital_period,
                        initial_eccentricity=initial_eccentricity,
                        initial_inclination=0.0,
                        disk_lock_frequency=self.Wdisk, 
                        disk_dissipation_age=self.disk_dissipation_age, 
                        secondary_formation_age=self.disk_dissipation_age) 
        
        binary.primary.select_interpolation_region(primary.core_formation_age())
        if isinstance(secondary, EvolvingStar):binary.secondary.detect_stellar_wind_saturation()

        binary.configure(age=primary.core_formation_age(),
                         semimajor=float('nan'),
                         eccentricity=float('nan'),
                         spin_angmom=numpy.array([0.0]),
                         inclination=None,
                         periapsis=None,
                         evolution_mode='LOCKED_SURFACE_SPIN')

        binary.primary.detect_stellar_wind_saturation()

        binary.secondary.configure(age=self.disk_dissipation_age, 
                            companion_mass=primary.mass,
                            semimajor=binary.semimajor(initial_orbital_period),
                            eccentricity=initial_eccentricity,
                            spin_angmom=spin_angmom,
                            inclination=inclination,
                            periapsis=periapsis,
                            locked_surface=False,
                            zero_outer_inclination=True,
                            zero_outer_periapsis=True
                            )


        return binary


    def initial_condition_errfunc(self,initial_conditions): 

        initial_orbital_period=initial_conditions[0]
        initial_eccentricity=initial_conditions[1]

        print('\nTrying Porb_initial = {} , e_initial = {}'.format(initial_orbital_period,initial_eccentricity))

        if initial_eccentricity>0.45 or initial_orbital_period<0 or initial_eccentricity<0:
            print('Invalid values')
            return scipy.nan,scipy.nan
       
        binary=self.create_binary_system(
            self.primary,
            self.secondary,
            initial_orbital_period,
            initial_eccentricity,
            secondary_angmom=self.secondary_angmom
            )
        
        binary.evolve(
            self.age, 
            self.evolution_max_time_step, 
            self.evolution_precision, 
            None,
            timeout=600
        )
                        
        final_state=binary.final_state()
        assert(final_state.age==self.age) 

        self.final_orbital_period=binary.orbital_period(final_state.semimajor)
        self.final_eccentricity=final_state.eccentricity

        self.delta_p=self.final_orbital_period-initial_orbital_period
        self.delta_e=self.final_eccentricity-initial_eccentricity

        if numpy.logical_or(numpy.isnan(self.delta_p),(numpy.isnan(self.delta_e))):
            print('Binary system was destroyed')
            raise ValueError

        self.spin=(2*numpy.pi*binary.primary.envelope_inertia(final_state.age)/final_state.primary_envelope_angmom)

        binary.delete()

        print('delta_p = {} , delta_e = {}'.format(self.delta_p,self.delta_e))
        print('Spin Period = ',self.spin)

        return self.delta_p,self.delta_e

    def initial_condition_solver(self):
        
        initial_guess=[self.Porb,self.eccentricity]
        #try:
        sol=scipy.optimize.root(self.initial_condition_errfunc,
                                initial_guess,
                                method='lm',
                                tol=1e-6,
                                options={'xtol':1e-6,
                                         'ftol':1e-6}

        )               

       # except Exception as e:
        #    print(e)
         #   print('Cannot calculate spin period')
         #   self.spin=scipy.nan
         #   return scipy.nan,scipy.nan

        initial_orbital_period_sol,initial_eccentricity_sol=sol.x

        print('Solver Results:')
        print('Intial Orbital Period = {} , Initial Eccentricity = {}'.format(initial_orbital_period_sol,initial_eccentricity_sol))
        print('Final Orbital Period = {} , Final Eccentricity = {}'.format(self.final_orbital_period,self.final_eccentricity))
        print('Errors: delta_p = {} , delta_e = {}'.format(self.delta_p,self.delta_e))
        print('Final Spin Period = {}'.format(self.spin))

        return initial_orbital_period_sol,initial_eccentricity_sol

    def initial_secondary_angmom(self):

        print('Calculating secondary angomom')
        star = self.create_star(self.secondary_mass,1)
        planet = self.create_planet(1.0)
        binary = self.create_binary_system(star,
                                      planet,
                                      10.0,
                                      0.0) 
        

        binary.evolve(self.disk_dissipation_age, 1e-3, 1e-6, None)

        disk_state = binary.final_state()

        planet.delete()
        star.delete()
        binary.delete()

        return numpy.array([disk_state.envelope_angmom, disk_state.core_angmom])


    def __init__(self,
                system,
                interpolator,
                sampling_parameters,
                fixed_parameters,
                mass_ratio,
                evolution_max_time_step=1e-3,
                evolution_precision=1e-6):

        self.system=system
        self.interpolator=interpolator
        self.sampling_parameters=sampling_parameters

        for item,value in sampling_parameters.items():
            setattr(self,item,value)

        for item,value in fixed_parameters.items():
            setattr(self,item,value)

        self.convective_phase_lag=phase_lag(self.logQ)
        self.mass_ratio=mass_ratio

        self.evolution_max_time_step=evolution_max_time_step
        self.evolution_precision=evolution_precision


    def __call__(self):


        self.secondary_mass=self.primary_mass*self.mass_ratio
        if numpy.logical_or((numpy.logical_or(self.secondary_mass>1.2,
                                              self.secondary_mass<0.4)),
                            (numpy.logical_or(numpy.isnan(self.primary_mass),
                                              numpy.isnan(self.secondary_mass)))
                            ):
            print('mass out of range')
            sys.stdout.flush()
            return scipy.nan

        self.secondary_angmom=self.initial_secondary_angmom()
        print('Secondary Spin Angmom = ',self.secondary_angmom)

        self.primary = self.create_star(self.primary_mass, 1)
        self.secondary = self.create_star(self.secondary_mass, 1)

        initial_orbital_period,initial_eccentricity=self.initial_condition_solver()

        self.primary.delete()
        self.secondary.delete()

        solver_results=dict()

        solver_results.update(self.sampling_parameters)
        solver_results['initial_orbital_period']=initial_orbital_period
        solver_results['initial_eccentricity']=initial_eccentricity
        solver_results['final_orbital_period']=self.final_orbital_period
        solver_results['final_eccentricity']=self.final_eccentricity
        solver_results['delta_p']=self.delta_p
        solver_results['delta_e']=self.delta_e
        solver_results['spin_period']=self.spin

        with open('solver_results_system_'+self.system+'.pickle','wb') as f:
            pickle.dump(solver_results,f)

        return self.spin





