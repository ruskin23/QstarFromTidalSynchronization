#!/usr/bin/env python3


import sys
sys.path.append('home/ruskin/projects/poet/PythonPackage')
sys.path.append('home/ruskin/projects/poet/scripts')

from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from orbital_evolution.binary import Binary
from orbital_evolution.transformations import phase_lag
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet
from initial_condition_solver import  InitialConditionSolver
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

    def create_star(self, mass):
        star = EvolvingStar(mass=mass,
                            metallicity=self.feh,
                            wind_strength=self.wind_strength if self.wind else 0.0,
                            wind_saturation_frequency=self.wind_saturation_frequency,
                            diff_rot_coupling_timescale=self.diff_rot_coupling_timescale,
                            interpolator=self.interpolator)

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

        for item,value in observational_parameters.items():
            setattr(self,item,value)


        self.convective_phase_lag=0.0

        self.target=Structure(age=self.age,
                              Porb=self.Porb,
                              Wdisk=self.Wdisk,
                              eccentricity=self.eccentricity)





    def initial_secodnary_angmom(self):

        star = self.create_star(self.secondary_mass)
        planet = self.create_planet(1.0)

        binary = self.create_binary_system(star,
                                      planet,
                                      10.0,
                                      self.disk_dissipation_age)

        binary.evolve(self.disk_dissipation_age, 1e-3, 1e-6, None)
        disk_state = binary.final_state()

        planet.delete()
        star.delete()
        binary.delete()

        print ('star-planet evolution completed')

        return numpy.array([disk_state.envelope_angmom, disk_state.core_angmom])

    def __call__(self,q,sol_file,frequency_break=None):

        self.convective_phase_lag=phase_lag(q)

        while True:
            IntialSecondaryAngmom=self.initial_secodnary_angmom()

            primary=self.create_star(self.primary_mass)
            secondary=self.create_star(self.secondary_mass)
            find_ic=InitialConditionSolver(system=self.system,
                                           print_cfile=self.print_cfile,
                                           breaks=self.breaks,
                                           disk_dissipation_age=self.disk_dissipation_age,
                                           evolution_max_time_step=1e-3,
                                           secondary_angmom=IntialSecondaryAngmom,
                                           is_secondary_star=True)


            print('Wdisk: ',self.Wdisk)
            print('secondary_angmom_initial = ',IntialSecondaryAngmom)
            solutions = find_ic(target=self.target, primary=primary,secondary=secondary)

            if bool(solutions)==False:
                self.Wdisk=numpy.random.uniform(low=numpy.pi/14,
                                                 high=numpy.pi/1.4)
                continue
            else:break

        print('Solution: ',solutions )

        sol=[]

        for key,value in solutions.items():
            sol.append(repr(value))

        sol='\t'.join(sol)

        if self.breakPower<0:q=logQ1

        with open(sol_file,'a',1) as f:
            f.write(repr(q)+'\t'+sol+'\n')

        primary.delete()
        secondary.delete()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('index',help='select system to run')
    parser.add_argument('-b',action='store',dest='breaks',
                        help='decide if breaks needed or not')
    parser.add_argument('-a',action = 'store_const',dest='add',const='add',help='add to logQvsPspin')
    parser.add_argument('-n',action = 'store_const',dest='new',const='new',help='make new table')
    args = parser.parse_args()

    serialized_dir ="/home/ruskin/projects/poet/stellar_evolution_interpolators"
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        b"eccentricity_expansion_coef.txt"
    )


    system=args.index
    print('System = ' ,system)

    data_file='SpinlogQCatalog_el0.4.txt'

    parameters=dict()

    spin_vs_logQ_file='break'+args.breaks+'/SpinLogQ_'+system+'.txt'

    if args.new:action='w'
    if args.add:action='a'

    breakPower=float(args.breaks)

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
            at_system=x[0]
            if at_system==system:
                parameters['primary_mass']=float(x[15])
                parameters['age']=float(x[16])
                parameters['feh']=float(x[17])

                parameters['eccentricity']=float(x[8])
                parameters['Porb']=float(x[6])
                mass_ratio=float(x[14])
                parameters['secondary_mass']=parameters['primary_mass']*mass_ratio
                parameters['Pspin']=float(x[12])


                parameters['Wdisk']=4.1
                parameters['disk_dissipation_age']=5e-3
                parameters['wind']=True
                parameters['wind_saturation_frequency']=2.54
                parameters['diff_rot_coupling_timescale']=5e-3
                parameters['wind_strength']=0.17

                parameters['system']=system
                parameters['print_cfile']=False
                parameters['breaks']=breakPower

                print('Mass Ratio = ', mass_ratio)

                logQ=numpy.linspace(5.0,6.5,5)
                #logQ=[8.2,8.3,8.4,8.5]

                print(breakPower)

                for q in logQ:

                    if breakPower==0.0:
                        TidalFrequencyBreaks=None
                        TidalFrequencyPowers=numpy.array([0.0])
                    elif breakPower>0.0:
                        TidalFrequencyBreaks=numpy.array([2*numpy.pi])
                        TidalFrequencyPowers=numpy.array([breakPower,breakPower])
                    elif breakPower<0.0:

                        #Calculate reference frequency
                        logQMax=4.0
                        PhaseLagMax=phase_lag(logQMax)
                        logQ1=q
                        PhaseLag1=phase_lag(logQ1)
                        omega1=2*numpy.pi
                        omegaRef=omega1*((PhaseLagMax/PhaseLag1)**(1.0/breakPower))

                        #set logQ to logQMax and define logQ1 parameter
                        parameters['logQ1']=q
                        q=logQMax

                        TidalFrequencyBreaks=numpy.array([omegaRef])
                        TidalFrequencyPowers=numpy.array([0.0,breakPower])

                    parameters['breakPower']=breakPower
                    parameters['tidal_frequency_breaks']=TidalFrequencyBreaks
                    parameters['tidal_frequency_powers']=TidalFrequencyPowers

                    print('parameters:', parameters)

                    print('\nCalculating for logQ = ', q)
                    evolve = evolution(interpolator,parameters)
                    evolve(q,spin_vs_logQ_file)

                break

