#!/usr/bin/env python3
import sys
import os
from pathlib import Path
from directories import directories
import numpy
home_dir=str(Path.home())
path=directories(home_dir)
sys.path.append(path.poet_path+'/PythonPackage')
sys.path.append(path.poet_path+'/scripts')

from orbital_evolution.transformations import phase_lag

from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as \
    orbital_evolution_library


from orbital_evolution.binary import Binary
from orbital_evolution.transformations import phase_lag
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet
from astropy import constants

import argparse

class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)

def create_planet(mass=(constants.M_jup / constants.M_sun).to('')):
    """Return a configured planet to use in the evolution."""

    planet = LockedPlanet(mass=mass, radius=(constants.R_jup / constants.R_sun).to(''))
    return planet

def create_star(mass,
                parameters,
                interpolator):

    star = EvolvingStar(mass=mass,
                        metallicity=parameters.feh,
                        wind_strength=parameters.wind_strength,
                        wind_saturation_frequency=parameters.wind_saturation_frequency,
                        diff_rot_coupling_timescale=parameters.diff_rot_coupling_timescale,
                        interpolator=interpolator)

    star.set_dissipation(zone_index=0,
                         tidal_frequency_breaks=parameters.tidal_frequency_breaks,
                         spin_frequency_breaks=parameters.spin_frequency_breaks,
                         tidal_frequency_powers=parameters.tidal_frequency_powers,
                         spin_frequency_powers=parameters.spin_frequency_powers,
                         reference_phase_lag=parameters.phase_lag_max)
        

    return star


def create_binary_system(primary,
                         secondary,
                         disk_dissipation_age=5e-3,
                         disk_lock_frequency=5e-3,
                         initial_orbital_period=10.0,
                         initial_eccentricity=0.0,
                         secondary_angmom=numpy.array([0.1,0.1])):
    """Create a binary system to evolve from the given objects."""


    if isinstance(secondary, LockedPlanet):
        spin_angmom = numpy.array([0.0])
        inclination = None
        periapsis = None
    else:
        secondary.select_interpolation_region(disk_dissipation_age)
        spin_angmom = secondary_angmom
        inclination = numpy.array([0.0])
        periapsis = numpy.array([0.0])

    primary.select_interpolation_region(primary.core_formation_age())

    binary = Binary(primary=primary,
                    secondary=secondary,
                    initial_orbital_period=initial_orbital_period,
                    initial_eccentricity=initial_eccentricity,
                    initial_inclination=0.0,
                    disk_lock_frequency=disk_lock_frequency,
                    disk_dissipation_age=disk_dissipation_age,
                    secondary_formation_age=disk_dissipation_age)
    binary.configure(age=primary.core_formation_age(),
                        semimajor=float('nan'),
                        eccentricity=float('nan'),
                        spin_angmom=numpy.array([0.0]),
                        inclination=None,
                        periapsis=None,
                        evolution_mode='LOCKED_SURFACE_SPIN')


    secondary.configure(age=disk_dissipation_age,
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

    primary.detect_stellar_wind_saturation()
    if isinstance(secondary, EvolvingStar):secondary.detect_stellar_wind_saturation()

    return binary


def get_binary_system(interpolator,
                      parameters,
                      secondary_angmom=None):

    primary=create_star(parameters.primary_mass,
                        parameters,
                        interpolator)
    secondary=create_star(parameters.secondary_mass,
                          parameters,
                          interpolator)
    
    
    binary_system=create_binary_system(primary,
                                       secondary,
                                       disk_dissipation_age=parameters.disk_dissipation_age,
                                       disk_lock_frequency=parameters.Wdisk,
                                       initial_orbital_period=parameters.orbital_period,
                                       initial_eccentricity=parameters.eccentricity,
                                       secondary_angmom=secondary_angmom)
    
    return binary_system



def evolve_binary(interpolator,
                  parameters):


    if type(parameters) is dict:
        parameters=Struct(**parameters)


    star = create_star(parameters.secondary_mass,
                       parameters,
                       interpolator)
    planet = create_planet(1.0)
    binary = create_binary_system(star,
                                  planet,
                                  disk_dissipation_age=parameters.disk_dissipation_age,
                                  disk_lock_frequency=parameters.Wdisk,
                                  initial_orbital_period=10.0,
                                  initial_eccentricity=0.0)

    binary.evolve(parameters.disk_dissipation_age, 1e-3, 1e-6, None)
    disk_state = binary.final_state()
    secondary_angmom=numpy.array([disk_state.envelope_angmom, disk_state.core_angmom])

    planet.delete()
    star.delete()
    binary.delete()

    binary = get_binary_system(interpolator, parameters, secondary_angmom=secondary_angmom)

    for dt in [1/10**k for k in [2,3,4]]:
        print('\nStarting evolution with dt {!r}'.format(dt))
        binary = get_binary_system(interpolator, parameters, secondary_angmom=secondary_angmom)
        binary.evolve(parameters.age,
                    dt,
                    parameters.evolution_precision,
                    None,
                    timeout=3600)

        final_state=binary.final_state()
        print('Reached age {!r}'.format(final_state.age))
        if final_state.age==parameters.age:break
        else:print('Evolution did not reach required age. Reducing by 1e-1')

    final_orbital_period=binary.orbital_period(final_state.semimajor)
    final_eccentricity=final_state.eccentricity
    print('final_eccntricity {} \nfinal orbital period {} \nfinal age {}'.format(final_eccentricity, final_orbital_period, final_state.age))

    return binary

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-l',
                        dest='logfile',
                        default=None,
                        help='select logfile to get parameters')
    args = parser.parse_args()
    
    parameters=dict()


    parameters['disk_dissipation_age']=5e-3
    parameters['wind_saturation_frequency']=2.54
    parameters['diff_rot_coupling_timescale']=5e-3
    parameters['wind_strength']=0.17
    parameters['evolution_max_time_step']=1e-2
    parameters['evolution_precision']=1e-6
    parameters['inclination']=0.0
    parameters['spin_frequency_breaks']=None
    parameters['spin_frequency_powers']=numpy.array([0.0])

    parameters["primary_mass"]=0.9241915970492054
    parameters["secondary_mass"]=0.8593663624245791
    parameters["feh"]=0.46458608151644215
    parameters["age"]=6.289791713929844
    parameters["Wdisk"]=1.1954168402393623
    parameters["phase_lag_max"]=6.251364607650211e-09
    parameters["tidal_frequency_breaks"]=numpy.array([0.12566371,7.02469755])
    parameters["tidal_frequency_powers"]=numpy.array([0.0,1.51411758, 0.0])

    parameters['orbital_period']= 11.653641112734459
    parameters['eccentricity']= 0.4406876762698404

    print('\nParameters = ')
    print(parameters)
    parameters=Struct(**parameters)

    serialized_dir = path.poet_path +  "/stellar_evolution_interpolators"
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    eccentricity_path=os.path.join(path.poet_path,'eccentricity_expansion_coef_O400.sqlite').encode('ascii')

    orbital_evolution_library.prepare_eccentricity_expansion(
        eccentricity_path,
        1e-4,
        True,
        True
    )

    evolved_binary=evolve_binary(interpolator,parameters)


