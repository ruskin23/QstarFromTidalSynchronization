#!/usr/bin/env python3 -u

import argparse
import scipy
import sys
import os
import os.path

from pathlib import Path
home_dir=str(Path.home())

git_dir='/QstarFromTidalSynchronization/MCMC/combined'

if home_dir=='/home/rxp163130':
    poet_path=home_dir+'/poet/'
    current_directory=home_dir+git_dir
    samples_directory=home_dir+'/QstarFromTidalSynchronization/MCMC/mcmc_mass_age/samples/updated_samples'

if home_dir=='/home/ruskin':
    poet_path=home_dir+'/projects/poet/'
    current_directory=home_dir+'/projects'+git_dir
    samples_directory=home_dir+'/projects/QstarFromTidalSynchronization/MCMC/mcmc_mass_age/samples/updated_samples'

if home_dir=='/home1/06850/rpatel23':
    work_dir='/work/06850/rpatel23/stampede2'
    poet_path=work_dir+'/poet'
    current_directory=work_dir+git_dir
    samples_directory=work_dir+'/QstarFromTidalSynchronization/MCMC/mcmc_mass_age/samples/updated_samples'

sys.path.append(poet_path+'PythonPackage')
sys.path.append(poet_path+'scripts')

from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from metropolis_hasting_test import MetropolisHastings


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

    parser.add_argument('-b',
                        action='store',
                        dest='breaks',
                        help='flag this to include tidal frequency breaks'
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

    if home_dir=='/home/rxp163130':output_directory=current_directory+'/ganymede/MCMC_'+system_number+'/'
    if home_dir=='/home/ruskin':output_directory=current_directory+'/MCMC_TESTS/TEST5/'
    if home_dir=='/home1/06850/rpatel23':output_directory=current_directory+'/stampede2/MCMC_'+system_number+'/'
    if os.path.isdir(output_directory)==False:os.mkdir(output_directory)

    catalog_file=current_directory+'/SpinlogQCatalog_el0.4.txt'
    solution_file=current_directory+'/SolutionFileBreaks0.0.txt'
    samples_file=samples_directory+'/MassAgeFehSamples_'+system_number+'.txt'
    stepfilename = output_directory+'step_file_'+args.instance+'.txt'

    with open(catalog_file,'r') as f:
        next(f)
        for lines in f:
            x=lines.split()
            at_system=x[0]
            if system_number==at_system:
                Porb_value=float(x[6])
                Porb_error=float(x[7])
                eccentricity_value=float(x[8])
                eccentricity_error=float(x[9])
                Pspin_value=float(x[12])
                Pspin_error=float(x[13])
                mass_ratio=float(x[14])
                primary_mass_value=float(x[15])
                age_value=float(x[16])
                feh_value=float(x[17])
                break


    with open(solution_file,'r') as f:
        next(f)
        for lines in f:
            x=lines.split()
            at_system=x[0]
            if system_number==at_system:
                logQ_value=float(x[1])
                break

    sampling_parameters=dict(Porb=dict(value=Porb_value,
                                       sigma=Porb_error,
                                       dist='Normal',
                                       step=Porb_error),

                             eccentricity=dict(value=eccentricity_value,
                                               sigma=eccentricity_error,
                                               dist='Normal',
                                               step=eccentricity_error),

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

                             primary_mass=dict(value=primary_mass_value,
                                               dist='Samples',
                                               step=1.0),

                             age=dict(value=age_value,
                                      dist='Samples',
                                      step=1.0),

                             feh=dict(value=feh_value,
                                      dist='Samples',
                                      step=0.3)

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


test_case='correlated_logQ'
mcmc = MetropolisHastings(system_number,
                          interpolator,
                          sampling_parameters,
                          fixed_parameters,
                          observed_spin,
                          catalog_file,
                          solution_file,
                          samples_file,
                          mass_ratio,
                          instance,
                          output_directory,
                          test_case)

sys.stdout.flush()
if args.start: mcmc.iterations()
elif args.cont: mcmc.continue_last()
else: print('provide correct arguments')

