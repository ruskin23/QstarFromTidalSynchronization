#!/usr/bin/env python3


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
from Solver1 import  InitialConditionSolver
from intial_secondary_angmom import IntialSecondaryAngmom
from basic_utils import Structure
import numpy
import scipy
from astropy import units, constants
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
                                tidal_frequency_breaks=self.tidal_frequency_breaks,
                                spin_frequency_breaks=None,
                                tidal_frequency_powers=self.tidal_frequency_powers,
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



    def __init__(self,interpolator,observational_parameters):

        self.interpolator=interpolator

        self.parameters=parameters

        for item,value in observational_parameters.items():
            setattr(self,item,value)


        self.convective_phase_lag=0.0

        self.target=Structure(age=self.age,
                              Porb=self.Porb,
                              Wdisk=self.Wdisk,
                              eccentricity=self.eccentricity)

    def __call__(self,q,sol_file,option=None):

        tdisk = self.disk_dissipation_age

        self.convective_phase_lag=phase_lag(q)

        secondary_angmom=IntialSecondaryAngmom(self.interpolator,
                                               q,
                                               self.parameters
                                               )


        primary = self.create_star(self.primary_mass, 1)
        secondary = self.create_star(self.secondary_mass, 1)
        find_ic = InitialConditionSolver(disk_dissipation_age=tdisk,
                                         evolution_max_time_step=1e-3,
                                         secondary_angmom=secondary_angmom,
                                         is_secondary_star=True)



        solutions = find_ic(target=self.target, primary=primary,secondary=secondary)

        if self.FindCircularLimit==True:
            return solutions['sol_e']

        if self.FindSyncLimit==True:
            return solution['spin']

        sol=[]

        for key,value in solutions.items():
            sol.append(repr(value))

        sol='\t'.join(sol)
        with open(sol_file,'a',1) as f:
            f.write(repr(q)+'\t'+sol+'\n')

        primary.delete()
        secondary.delete()



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('system',help='select system to run')
    parser.add_argument('-a',action = 'store_const',dest='add',const='add',help='add to logQvsPspin')
    parser.add_argument('-n',action = 'store_const',dest='new',const='new',help='make new table')
    parser.add_argument('-b',action = 'store_const',dest='breaks',const='breaks',help='breaks')
    args = parser.parse_args()

    serialized_dir ="/home/ruskin/projects/poet/stellar_evolution_interpolators"
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        b"eccentricity_expansion_coef.txt"
    )

    system=int(args.system)

    print('System = ' ,system)

    data_file='spin_vs_logQ_systems.txt'
    parameters=dict()

    if args.breaks:spin_vs_logQ_file='../break2.0/SpinLogQ_WithBreaks_'+args.system+'.txt'
    else: spin_vs_logQ_file='../SpinLogQFiles/SpinLogQ_'+args.system+'.txt'

    if args.new:action='w'
    if args.add:action='a'

    with open(spin_vs_logQ_file,action) as f:
        f.write('logQ'+'\t'+
                'spin'+'\t'+
                'Porb_initial'+'\t'+
                'e_initial'+'\t'+
                'Porb_current'+'\t'+
                'e_current'+'\t'+
                'delta_p'+'\t'+
                'detla_e'+'\n')

    with open(data_file,'r') as f:
        next(f)
        for lines in f:
            x=lines.split()
            at_system=int(x[0])
            print(at_system)
            if at_system==system:
                parameters['primary_mass']=float(x[15])
                parameters['age']=10**(float(x[16]))
                parameters['feh']=float(x[17])

                parameters['eccentricity']=float(x[8])
                parameters['Porb']=float(x[6])
                parameters['Pspin']=float(x[12])
                mass_ratio=float(x[14])
                parameters['secondary_mass']=parameters['primary_mass']*mass_ratio

                if args.breaks:
                    TidalFrequencyBreaks=numpy.array([abs(-4*numpy.pi*((1.0/parameters['Porb'])
                                                                  -
                                                                  1.0/parameters['Pspin']))])
                    TidalFrequencyPowers=numpy.array([2.0,2.0])
                    parameters['breaks']=True
                else:
                    TidalFrequencyBreaks=None
                    TidalFrequencyPowers=numpy.array([0.0])
                    parameters['breaks']=False

                parameters['tidal_frequency_breaks']=TidalFrequencyBreaks
                parameters['tidal_frequency_powers']=TidalFrequencyPowers

                parameters['Wdisk']=4.3
                parameters['disk_dissipation_age']=5e-3
                parameters['wind']=True
                parameters['wind_saturation_frequency']=2.54
                parameters['diff_rot_coupling_timescale']=5e-3
                parameters['wind_strength']=0.17

                parameters['FindCircularLimit']=False
                parameters['FindSyncLimit']=False
                print('Parameters: ', parameters)

                evolve = evolution(interpolator,parameters)

                logQ = numpy.linspace(9.1,10.0,5)
                #logQ=[5.775]
                for q in logQ:
                    print('\n\nCalculating for logQ = ', q)
                    evolve(q,spin_vs_logQ_file,option=1)

                break
