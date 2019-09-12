#!/usr/bin/env python3

import sys

from pathlib import Path
home_dir=str(Path.home())

if home_dir=='/home/rxp163130':
    poet_path=home_dir+'/poet'
if home_dir=='/home/ruskin':
    poet_path=home_dir+'/projects/poet/'

sys.path.append(poet_path+'PythonPackage')
sys.path.append(poet_path+'scripts')

import sys

from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from orbital_evolution.binary import Binary
from orbital_evolution.transformations import phase_lag
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet
from basic_utils import Structure
import numpy
import scipy
from astropy import units, constants
import pickle


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



        return star

    def create_binary_system(self,
                             primary,
                             secondary,
                             secondary_angmom=None):
        """Create a binary system to evolve from the given objects."""

        if isinstance(secondary, LockedPlanet):
            secondary_config = dict(spin_angmom=numpy.array([0.0]),
                                    inclination=None,
                                    periapsis=None)
        else:
            secondary.select_interpolation_region(self.disk_dissipation_age)
            secondary_config = dict(spin_angmom=secondary_angmom,
                                    inclination=numpy.array([0.0]),
                                    periapsis=numpy.array([0.0]))


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
                        secondary_formation_age=30.0)


        secondary.configure(age=30.0,
                            companion_mass=primary.mass,
                            semimajor=binary.semimajor(self.Porb),
                            eccentricity=0.0,
                            locked_surface=False,
                            zero_outer_inclination=True,
                            zero_outer_periapsis=True,
                            **secondary_config)

        binary.configure(age=primary.core_formation_age(),
                         semimajor=float('nan'),
                         eccentricity=float('nan'),
                         spin_angmom=numpy.array([0.0]),
                         inclination=None,
                         periapsis=None,
                         evolution_mode='LOCKED_SURFACE_SPIN')

        return binary

    def __init__(self,
                 interpolator,
                 fixed_parameters,
                 age,
                 primary_mass,
                 secondary_mass
                 ):

        self.interpolator=interpolator

        for item,value in fixed_parameters.items():
            setattr(self,item,value)

        self.age=age
        self.primary_mass=primary_mass
        self.secondary_mass=secondary_mass


    def __call__(self):


        star = self.create_star(self.secondary_mass)
        planet = self.create_planet(1.0)

        binary = self.create_binary_system(star,
                                      planet
                                           )

        binary.evolve(self.disk_dissipation_age, 1e-3, 1e-6, None)

        disk_state = binary.final_state()
        secondary_angmom=numpy.array([disk_state.envelope_angmom,disk_state.core_angmom])

        planet.delete()
        star.delete()
        binary.delete()

        print ('star-planet evolution completed')


        primary = self.create_star(self.primary_mass)
        secondary = self.create_star(self.secondary_mass)


        binary = self.create_binary_system(primary,
                                           secondary,
                                           secondary_angmom=secondary_angmom
                                           )
        binary.evolve(self.age,1e-3,1e-6,None)

        final_state = binary.final_state()
        assert(final_state.age==self.age)

        spin =  (2.0 * scipy.pi *
                binary.primary.envelope_inertia(final_state.age)
                /
                final_state.primary_envelope_angmom
        )

        primary.delete()
        secondary.delete()

        return spin


if __name__ == '__main__':

    serialized_dir = poet_path +  "stellar_evolution_interpolators"
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        b"eccentricity_expansion_coef.txt"
    )

    data_file='/mnt/md0/ruskin/QstarFromTidalSynchronization/binary_star_evolution/mass_calculations/teff_logg_feh.txt'
    spin_output_file='/mnt/md0/ruskin/QstarFromTidalSynchronization/binary_star_evolution/mass_calculations/spin_output.txt'

    with open(spin_output_file,'w') as f:
        f.write('###\n')

    with open(data_file,'r') as f:
        next(f)
        for k,line in enumerate(f):

            #if k>1:break
            print('On line ', k)

            data=line.split()

            if len(data)==3:
                print('no solution for this')
                continue

            star_parameters={'age' : [float(data[3])],
                        'primary_mass' : [float(data[4])],
                        'secondary_mass' : [float(data[5])]
                        }

            if len(data)==9:
                star_parameters['age'].append(float(data[6]))
                star_parameters['primary_mass'].append(float(data[7]))
                star_parameters['secondary_mass'].append(float(data[8]))
                two_solution=True
            print('Star parameters: ',star_parameters)

            teff=data[0]
            logg=data[1]
            feh=float(data[2])


            fixed_parameters = dict(
                                Porb=1.798,
                                Wdisk=4.48,
                                feh=feh,
                                disk_dissipation_age=5e-3,
                                wind=True,
                                wind_saturation_frequency=2.54,
                                diff_rot_coupling_timescale=5e-3,
                                wind_strength=0.17,
                                inclination=scipy.pi/2
                            )

            with open(spin_output_file,'a',1) as f:
                print('opened_file')

                for i in range(len(star_parameters['age'])):
                    print('Setting up evolution for age',star_parameters['age'][i])

                    get_evolution=evolution(interpolator,
                                        fixed_parameters,
                                        star_parameters['age'][i],
                                        star_parameters['primary_mass'][i],
                                        star_parameters['secondary_mass'][i]
                                        )
                    spin=get_evolution()
                    print('Evolution Complete')

                    f.write(teff+'\t'+
                        logg+'\t'+
                        repr(feh)+'\t'+
                        repr(star_parameters['age'][i])+'\t'+
                        repr(star_parameters['primary_mass'][i])+'\t'+
                        repr(star_parameters['secondary_mass'][i])+'\t'+
                        repr(spin)+'\n'
                        )



