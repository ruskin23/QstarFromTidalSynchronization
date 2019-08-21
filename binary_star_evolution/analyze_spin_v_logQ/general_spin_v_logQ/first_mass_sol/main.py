#!/usr/bin/env python3
import os

import sys
sys.path.append('.../poet/PythonPackage')
sys.path.append('.../poet/scripts')


from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from orbital_evolution.binary import Binary
from orbital_evolution.transformations import phase_lag
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet
from mass_calculation import Derive_mass
from initial_condition_solver import  InitialConditionSolver
from basic_utils import Structure
import numpy
import scipy
from astropy import units, constants
import pickle
import argparse
from scipy.optimize import brentq

class evolution:

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

        secondary.configure(age=self.disk_dissipation_age,
                            companion_mass=primary.mass,
                            semimajor=initial_semimajor,
                            eccentricity=0.0,
                            locked_surface=False,
                            zero_outer_inclination=True,
                            zero_outer_periapsis=True,
                            **secondary_config)

        print ("BEGINSAT")

        if isinstance(secondary, EvolvingStar):
            secondary.detect_stellar_wind_saturation()
            print ("DETECTED SECONDARY WIND SAT")



        primary.select_interpolation_region(primary.core_formation_age())
        primary.detect_stellar_wind_saturation()

        binary = Binary(primary=primary,
                        secondary=secondary,
                        initial_orbital_period=self.Porb,
                        initial_eccentricity=0.0,
                        initial_inclination=0.0,
                        disk_lock_frequency=self.Wdisk,
                        disk_dissipation_age=self.disk_dissipation_age,
                        secondary_formation_age=self.disk_dissipation_age)

        binary.configure(age=primary.core_formation_age(),
                         semimajor=float('nan'),
                         eccentricity=float('nan'),
                         spin_angmom=numpy.array([0.0]),
                         inclination=None,
                         periapsis=None,
                         evolution_mode='LOCKED_SURFACE_SPIN')

        return binary


    def calculate_star_masses(self):
        star_masses = []

        print ("Calculating Masses\n")

        mass = Derive_mass(
                            self.interpolator,
                            self.teff_primary,
                            self.logg,
                            self.feh)

        sol = mass()
        self.primary_mass = sol[0]
        self.age=sol[1]
        self.secondary_mass = self.primary_mass*self.mass_ratio

        print(self.primary_mass)
        print(self.secondary_mass)
        print(self.age)


    def __init__(self,interpolator,observational_parameters):

        self.interpolator=interpolator

        for item,value in observational_parameters.items():
            setattr(self,item,value)


        self.primary_mass=0.0
        self.secondary_mass=0.0
        self.calculate_star_masses()

        self.convective_phase_lag=0.0

        self.target=Structure(age=self.age,
                              Porb=self.Porb,
                              Wdisk=self.Wdisk,
                              eccentricity=self.eccentricity)

    def __call__(self,q,sol_file,option=None):

        tdisk = self.disk_dissipation_age

        self.convective_phase_lag=phase_lag(q)

        star = self.create_star(self.secondary_mass,1)
        planet = self.create_planet(1.0)

        binary = self.create_binary_system(star,
                                      planet,
                                      10.0,
                                      tdisk)

        binary.evolve(tdisk, 1e-3, 1e-6, None)

        disk_state = binary.final_state()


        planet.delete()
        star.delete()
        binary.delete()

        print ('star-planet evolution completed')

        primary = self.create_star(self.primary_mass, 1)
        secondary = self.create_star(self.secondary_mass, 1)
        find_ic = InitialConditionSolver(disk_dissipation_age=tdisk,
                                         evolution_max_time_step=1e-3,
                                         secondary_angmom=numpy.array(
                                             [disk_state.envelope_angmom, disk_state.core_angmom]),
                                         is_secondary_star=True)



        initial_p, initial_e, current_p, current_e, current_spin, delta_p, delta_e = find_ic(target=self.target, primary=primary,secondary=secondary)


        if option == 1:
            with open(sol_file,'a') as f:
                f.write(repr(q) + '\t' + repr(current_spin) + '\t' + repr(initial_p) + '\t'
                    + repr(initial_e) + '\t' + repr(current_p) + '\t' +
                    repr(current_e) + '\t' + repr(self.convective_phase_lag) +
                    '\t' + repr(delta_p) + '\t' + repr(delta_e) + '\n')



        return current_spin

        primary.delete()
        secondary.delete()

if __name__ == '__main__':

    #parser = argparse.ArgumentParser()
    #parser.add_argument('instance',help='select system to run')
    #args = parser.parse_args()

    serialized_dir ="/home/ruskin/projects/poet/stellar_evolution_interpolators"
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        b"eccentricity_expansion_coef.txt"
    )
    fsol_range=open('solution_range_file','w')
    fsol_range.write('KIC'+'\t'+'logq_min'+'\t'+'logq_max'+'\n')

    #data_file = 'catalog_'+args.instance+'.txt'
    files = [1,2,3,4,5,6,7,8,9,10]
    #files = [11,12,13,14,15,16,17,18,19,20]
    #files = [21,22,23,24,25,26,27,28,29,30]
    #files = [31,32,33,34,35,36,37]
    for i in files:
        data_file='catalog_'+repr(i)+'_p.txt'

        with open(data_file,'r') as f:
            #next(f)
            for lines in f:
                data=lines.split()
                KIC='KIC'+data[0]
                teff=float(data[1])
                feh=float(data[2])
                logg=float(data[3])
                eccentricity=float(data[4])
                Porb=float(data[5])
                Pspin=float(data[6])
                mass_ratio=float(data[7])

                print('printing parameters:')
                print('KIC: ',KIC)
                print('teff: ',teff)
                print('feh: ',feh)
                print('logg: ',logg)
                print('eccentricity: ',eccentricity)
                print('Porb: ',Porb)
                print('Pspin: ',Pspin)
                print('mass_ratio: ',mass_ratio)

                spin_vs_logQ_file='spin_vs_logQ_'+KIC+'.txt'
                if os.path.isfile(spin_vs_logQ_file)==True:
                    break
                with open(spin_vs_logQ_file,'w') as f_svq:
                    f_svq.write('#DATA:' + '\n' + '#KIC' + '\t' + 'Teff'  + '\t' +  'FeH' + '\t' + 'logg' + '\t' + 'eccentricity' + '\t' + 'Porb' + '\t' + 'Pspin' + '\t' + 'q' + '\n' + '#')
                    f_svq.write(lines)

                parameters = dict(
                        teff_primary=teff,
                        feh= feh,
                        logg=logg,
                        Wdisk=4.306699756301906,
                        Porb=Porb,
                        incination=0.0,
                        disk_dissipation_age=5e-3,
                        wind=True,
                        planet_formation_age=5e-3,
                        wind_saturation_frequency=2.54,
                        diff_rot_coupling_timescale=5e-3,
                        wind_strength=0.17,
                        eccentricity=eccentricity,
                        mass_ratio=mass_ratio)

                evolve = evolution(interpolator,parameters)

                check_sign=1
                q_max=0.0
                q_min=0.0
                logQ = numpy.arange(6.0,11.0,1.0)
                for q in logQ:
                    print('Calculating for logQ = ', q)
                    spin = evolve(q,spin_vs_logQ_file,option=1)
                    print('Obtained spin = ', spin)
                    spin_diff = Pspin-spin
                    print('spin difference = ',spin_diff)
                    check_sign=check_sign*spin_diff
                    if check_sign>0 and q_max==0:q_min=q
                    if check_sign<0 and q_max==0:q_max=q
                fsol_range.write(KIC+'\t'+repr(q_min)+'\t'+repr(q_max)+'\n')

#            logQ_sol = brentq(
#               lambda q:evolve(q,None,option=2)-Pspin,
#               q_max,
#               q_min)


#            with open('logQ_solution','w') as f:
#                f.write(repr(logQ_sol))




