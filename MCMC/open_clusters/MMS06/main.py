#!/usr/bin/env python3 -u

import argparse
import scipy
import sys
import os
import os.path

from pathlib import Path
home_dir=str(Path.home())

git_dir='/QstarFromTidalSynchronization/MCMC/open_clusters'

if home_dir=='/home/rxp163130':
    poet_path=home_dir+'/poet/'
    current_directory=home_dir+git_dir

if home_dir=='/home/ruskin':
    poet_path=home_dir+'/projects/poet/'
    current_directory=home_dir+'/projects'+git_dir

if home_dir=='/home1/06850/rpatel23':
    work_dir='/work/06850/rpatel23/stampede2'
    poet_path=work_dir+'/poet'
    current_directory=work_dir+git_dir

sys.path.append(poet_path+'PythonPackage')
sys.path.append(poet_path+'scripts')

import random

from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from metropolis_hasting import MetropolisHastings


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

    serialized_dir = poet_path +  "/stellar_evolution_interpolators"
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    eccentricity_path=os.path.join(poet_path,'eccentricity_expansion_coef.txt').encode('ascii')

    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        eccentricity_path
    )


    args = cmdline_args()
    instance = args.instance
    system_number=args.system

    if home_dir=='/home/rxp163130':output_directory=current_directory+'/Cluster_Results/ganymede/MCMC_'+system_number+'/'
    if home_dir=='/home/ruskin':output_directory=current_directory+'/kartof/MCMC_'+system_number+'/'
    if home_dir=='/home1/06850/rpatel23':output_directory=current_directory+'/Cluster_Results/stampede2/MCMC_'+system_number+'/'

    catalog_file=current_directory+'/catalog.txt'
    stepfilename = output_directory+'step_file_'+args.instance+'.txt'

    with open(catalog_file,'r') as f:
        next(f)
        for lines in f:
            x=lines.split()
            at_system=x[0]
            if system_number==at_system:
                Porb_value=float(x[1])
                Porb_error=1e-5
                eccentricity_value=float(x[2])
                eccentricity_error=float(x[3])
                primary_mass=float(x[4])
                primary_mass_error=0.1
                mass_functino=float(x[5])
                Pspin_value=float(x[6])
                Pspin_error=float(x[7])
                break

    age_value=0.1
    age_error=0.01
    feh_value=-0.21
    feh_error=0.1

    sampling_parameters=dict(Porb=dict(value=Porb_value,
                                       sigma=Porb_error,
                                       dist='Normal',
                                       step=Porb_error),

                             eccentricity=dict(value=eccentricity_value,
                                               sigma=eccentricity_error,
                                               dist='Normal',
                                               step=eccentricity_error),

                             primary_mass=dict(value=primary_mass_value,
                                               sigma=primary_mass_error,
                                               dist='Normal',
                                               step=1.0),

                             age=dict(value=age_value,
                                      sigma=age_error,
                                      dist='Normal',
                                      step=),

                             feh=dict(value=age_value,
                                      sigma=age_error,
                                      dist='Nomal',
                                      step=),

                             cosi=dict(value=cosi_value,
                                       sigma=cosi_error,
                                       dist='Normal',
                                       step=0.1),

                             Wdisk=dict(value=4.1,
                                        min=2*scipy.pi/14,
                                        max=2*scipy.pi/1.4,
                                        dist='Uniform',
                                        step=0.5),

                             logQ=dict(value=logQ_value,
                                       min=5.0,
                                       max=12.0,
                                       dist='Uniform',
                                       step=0.5),

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


with open(stepfilename,'w') as f:
    for key,value in sampling_parameters.items():
        f.write(key + '\t' + repr(value['step']) + '\n')


mcmc = MetropolisHastings(system_number,
                          interpolator,
                          sampling_parameters,
                          fixed_parameters,
                          mass_function,
                          observed_spin,
                          instance,
                          current_directory,
                          output_directory)

sys.stdout.flush()
if args.start: mcmc.iterations()
elif args.cont: mcmc.continue_last()
else: print('provide correct arguments')

