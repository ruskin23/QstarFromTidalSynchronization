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

    print('\ncreated star with: \nmass = {} \nmetalliciy = {} \nwind_srength = {} \nwind_saturation_frequency = {} \ndiff_rot_coupling_timescale = {}'.format(
            repr(mass),
            repr(parameters.feh),
            repr(parameters.wind_strength),
            repr(parameters.wind_saturation_frequency),
            repr(parameters.diff_rot_coupling_timescale))
            )

    # star.select_interpolation_region(star.core_formation_age())

    star.set_dissipation(zone_index=0,
                         tidal_frequency_breaks=parameters.tidal_frequency_breaks,
                         spin_frequency_breaks=parameters.spin_frequency_breaks,
                         tidal_frequency_powers=parameters.tidal_frequency_powers,
                         spin_frequency_powers=parameters.spin_frequency_powers,
                         reference_phase_lag=parameters.phase_lag_max)
        
    print('\nset dissipation with \ntidal_frequency_breaks = {} \nspin_frequency_breaks = {} \ntidal_frequency_powers = {} \nspin_frequency_powers = {} \nreference_phase_lag = {}'.format(
            repr(parameters.tidal_frequency_breaks),
            repr(parameters.spin_frequency_breaks),
            repr(parameters.tidal_frequency_powers),
            repr(parameters.spin_frequency_powers),
            repr(parameters.phase_lag_max)))

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
        print('secondary angmom {}'.format(secondary_angmom))
        inclination = numpy.array([0.0])
        periapsis = numpy.array([0.0])

    primary.select_interpolation_region(primary.core_formation_age())

    print('\nCreating binary with:')
    print('\ninitial_orbital_period = {}'.format(repr(initial_orbital_period)))
    print('\ninitial_eccentricity = {}'.format(repr(initial_eccentricity)))
    print('\ndisk_lock_frequency = {}'.format(repr(disk_lock_frequency)))
    print('\ndisk_dissipation_age = {}'.format(repr(disk_dissipation_age)))
    print('\nseconday_formation_age = {}'.format(repr(disk_dissipation_age)))
    binary = Binary(primary=primary,
                    secondary=secondary,
                    initial_orbital_period=initial_orbital_period,
                    initial_eccentricity=initial_eccentricity,
                    initial_inclination=0.0,
                    disk_lock_frequency=disk_lock_frequency,
                    disk_dissipation_age=disk_dissipation_age,
                    secondary_formation_age=disk_dissipation_age)

    print('\nconfiguring binary with:')
    print('age = {}'.format(repr(primary.core_formation_age())))
    binary.configure(age=primary.core_formation_age(),
                        semimajor=float('nan'),
                        eccentricity=float('nan'),
                        spin_angmom=numpy.array([0.0]),
                        inclination=None,
                        periapsis=None,
                        evolution_mode='LOCKED_SURFACE_SPIN')

    print('\nconfiguring secondary with:')
    print('\nage = {}'.format(repr(disk_dissipation_age)))
    print('\ncompanion_mass = {}'.format(repr(primary.mass)))
    print('\nsemimajor = {}'.format(repr(binary.semimajor(initial_orbital_period))))
    print('\spin_angmom = {}'.format(repr(spin_angmom)))

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

    print('\nCreating Primary Star')
    primary=create_star(parameters.primary_mass,
                        parameters,
                        interpolator)
    print('\nCreating Secondary Star')
    secondary=create_star(parameters.secondary_mass,
                          parameters,
                          interpolator)
    
    print('\nevolving binary with \ndisk_dissipation_age = {} \ndisk_lock_frequency = {} \ninitial_orbital_period = {} \neccentricity = {}'.format(
           repr(parameters.disk_dissipation_age),
           repr(parameters.Wdisk),
           repr(parameters.orbital_period),
           repr(parameters.eccentricity)
    ))

    print('\nprimary_mass = {} \nprimary_feh = {} \nsecondary_mass = {} \nsecondary_feh = {}'.format(
        repr(primary.mass),
        repr(primary.metallicity),
        repr(secondary.mass),
        repr(secondary.metallicity)
    ))

    print('\nprimary_core_formation_age = {} \nsecondary_core_formation_age = {}'.format(
        repr(primary.core_formation_age()),
        repr(secondary.core_formation_age())
    ))
    
    binary_system=create_binary_system(primary,
                                       secondary,
                                       disk_dissipation_age=parameters.disk_dissipation_age,
                                       disk_lock_frequency=parameters.Wdisk,
                                       initial_orbital_period=parameters.orbital_period,
                                       initial_eccentricity=parameters.eccentricity,
                                       secondary_angmom=secondary_angmom)
    
    return binary_system



def evolve_binary(interp,
                  parameters):


    eccentricity_path=os.path.join(path.poet_path,'eccentricity_expansion_coef_O400.sqlite').encode('ascii')

    orbital_evolution_library.prepare_eccentricity_expansion(
        eccentricity_path,
        1e-4,
        True,
        True
    )

    serialized_dir = path.poet_path +  "/stellar_evolution_interpolators"
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')


    if type(parameters) is dict:
        print('\n Converting Parameters Dictionary to Struct')
        print(parameters)
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

    binary.evolve(parameters.age,
                  parameters.evolution_max_time_step,
                  parameters.evolution_precision,
                  None,
                  timeout=3600)

    final_state=binary.final_state()

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
    parameters['evolution_max_time_step']=1e-3
    parameters['evolution_precision']=1e-6
    parameters['inclination']=0.0
    parameters['spin_frequency_breaks']=None
    parameters['spin_frequency_powers']=numpy.array([0.0])

    parameters['Wdisk']=3.198953282185807
    parameters['primary_mass']=0.8220187083597497
    parameters['secondary_mass']=0.5384657839157082
    parameters['feh']=-0.5054652098816972
    parameters['age']=10.957202490925933

    parameters['phase_lag_max']=2.034701463188306e-07
    parameters['tidal_frequency_breaks']=numpy.array([0.12566371,0.6512592])
    parameters['tidal_frequency_powers']=numpy.array([0.0,1.18840619, 0.0])

    parameters['orbital_period']=15.233962642170873
    parameters['eccentricity']=0.7

    print(parameters)
    parameters=Struct(**parameters)
    
    serialized_dir = path.poet_path +  "/stellar_evolution_interpolators"
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    # eccentricity_path=os.path.join(path.poet_path,'eccentricity_expansion_coef_O400.sqlite').encode('ascii')

    # orbital_evolution_library.prepare_eccentricity_expansion(
    #     eccentricity_path,
    #     1e-4,
    #     True,
    #     True
    # )

    evolved_binary=evolve_binary(interpolator,parameters)



    
    # if args.logfile is not None:
    #     _simple_quantities=['primary_mass',
    #                 'secondary_mass',
    #                 'feh',
    #                 'age',
    #                 'Wdisk',
    #                 'orbital_period',
    #                 'eccentricity',
    #                 'phase_lag_max']


    #     with open('logfile/files/'+args.logfile,'r') as f:
    #         for i,lines in enumerate(f):

    #             if i>6 and i<42 and i%2==1:

    #                 x=lines.split()

    #                 if x[0] in _simple_quantities:
    #                     parameters[x[0]]=float(x[1])

    #                 if x[0]=='tidal_frequency_breaks':
    #                     x=lines.split()
    #                     if len(x)==2:
    #                         parameters[x[0]]=numpy.atleast_1d(float(x[1][1:-1]))
    #                     elif len(x)==3:
    #                         a=float(x[1][1:])
    #                         b=float(x[-1][:-1])
    #                         parameters[x[0]]=numpy.array([a,b])
    #                     elif len(x)==4:
    #                         a=float(x[1][1:])
    #                         b=float(x[-2])
    #                         parameters[x[0]]=numpy.array([a,b])

    #                 if x[0]=='tidal_frequency_powers':
    #                     if len(x)==4:
    #                         a=float(x[2])
    #                         b=float(x[3][0:-1])
    #                         value=numpy.array([a,b])
    #                     elif len(x)==5:
    #                         a=float(x[1][1:])
    #                         b=float(x[2])
    #                         c=float(x[3])
    #                         value=numpy.array([a,b,c])
    #                     parameters[x[0]]=numpy.array(value)
