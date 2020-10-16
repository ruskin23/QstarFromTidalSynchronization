#!/usr/bin/env python3 -u

import sys
import os

from pathlib import Path
from directories import directories

home_dir=str(Path.home())
path=directories(home_dir)
sys.path.append(path.poet_path+'/PythonPackage')
sys.path.append(path.poet_path+'/scripts')

from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from metropolis_hasting import MetropolisHastings

import argparse
import scipy


def cmdline_args():

    parser = argparse.ArgumentParser()

    parser.add_argument('-s',
                        action='store_const',
                        dest='start',
                        const='start',
                        help='start mcmc from beginning'
                        )

    parser.add_argument('-c',
                        action='store_const',
                        dest='cont',
                        const='continue',
                        help='continue mcmc from last iteration'
                        )

    parser.add_argument('-i',
                        action='store',
                        dest='instance',
                        help='define an instance of mcmc'
                        )

    parser.add_argument('-l',
                        action='store',
                        dest='system',
                        help='select a system for mcmc'
                        )


    return parser.parse_args()



if __name__ == '__main__':

    serialized_dir = path.poet_path +  "/stellar_evolution_interpolators"
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    eccentricity_path=os.path.join(path.poet_path,'eccentricity_expansion_coef.txt').encode('ascii')

    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        eccentricity_path
    )


    args = cmdline_args()
    instance = args.instance
    system_number=args.system

    output_directory=path.output_directory
    
    catalog_file=path.current_directory+'/WindeCatalog.txt'

    with open(catalog_file,'r') as f:
        next(f)
        for lines in f:
            x=lines.split()
            at_system=x[0]
            if system_number==at_system:
                Porb_value=float(x[2])
                Porb_error=float(x[3])
                ecosw_value=float(x[4])
                ecosw_error=float(x[5])
                esinw_value=float(x[6])
                esinw_error=float(x[7])
                primary_mass_value=float(x[8])
                primary_mass_sigma=float(x[9])
                secondary_mass_value=float(x[10])
                secondary_mass_sigma=float(x[11])
                age_value=float(x[12])
                age_sigma=float(x[13])
                feh_value=float(x[14])
                feh_sigma=float(x[15])
                Pspin_value=float(x[16])
                Pspin_error=float(x[17])
                break


    sampling_parameters=dict(Porb=dict(value=Porb_value,
                                       sigma=Porb_error,
                                       min_value=0.0,
                                       max_value=100.0,
                                       dist='Normal',
                                       step=Porb_error),

                             ecosw=dict(value=ecosw_value,
                                        sigma=ecosw_error,
                                        min_value=-0.5,
                                        max_value=0.5,
                                        dist='Normal',
                                        step=2*ecosw_error),

                             esinw=dict(value=esinw_value,
                                        sigma=esinw_error,
                                        min_value=-0.5,
                                        max_value=0.5,
                                        dist='Normal',
                                        step=2*esinw_error),


                             primary_mass=dict(value=primary_mass_value,
                                               sigma=primary_mass_sigma,
                                               min_value=0.4,
                                               max_value=1.2,
                                               dist='Normal',
                                               step=2*primary_mass_sigma),

                             secondary_mass=dict(value=secondary_mass_value,
                                                 sigma=secondary_mass_sigma,
                                                 min_value=0.4,
                                                 max_value=1.2,
                                                 dist='Normal',
                                                 step=2*secondary_mass_sigma),

                             age=dict(value=age_value,
                                      sigma=age_sigma,
                                      min_value=1e-3,
                                      max_value=12.0,
                                      dist='Normal',
                                      step=2*age_sigma
                                      ),

                             feh=dict(value=feh_value,
                                      sigma=feh_sigma,
                                      min_value=-1.014,
                                      max_value=0.537,
                                      dist='Normal',
                                      step=2*feh_sigma),

                             Wdisk=dict(value=2.75,
                                        min_value=2*scipy.pi/14,
                                        max_value=2*scipy.pi/1.4,
                                        dist='Uniform',
                                        step=1.0),

                             logQ=dict(value=7.0,
                                       min_value=5.0,
                                       max_value=12.0,
                                       dist='Uniform',
                                       step=1.0),

                             )


    observed_spin=dict(value=Pspin_value,
                       sigma=Pspin_error)

    fixed_parameters = dict(disk_dissipation_age=5e-3,
                            planet_formation_age=5e-3,
                            wind=True,
                            wind_saturation_frequency=2.54,
                            diff_rot_coupling_timescale=5e-3,
                            wind_strength=0.17,
                            inclination=scipy.pi/2
                            )

    print('Sampling Parameters: ',sampling_parameters)
    print('Observed Spin: ',observed_spin)



mcmc = MetropolisHastings(system_number,
                          interpolator,
                          sampling_parameters,
                          fixed_parameters,
                          observed_spin,
                          instance,
                          path.current_directory,
                          path.output_directory)

sys.stdout.flush()
if args.start: mcmc.iterations()
elif args.cont: mcmc.continue_last()
else: print('provide correct arguments')

