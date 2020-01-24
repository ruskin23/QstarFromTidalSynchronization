#!/usr/bin/env python3 -u

import argparse
import scipy
import sys
import os
import os.path

from pathlib import Path
home_dir=str(Path.home())

git_dir='/QstarFromTidalSynchronization/MCMC/general'

if home_dir=='/home/rxp163130':
    poet_path=home_dir+'/poet/'
    current_directory=home_dir+git_dir

if home_dir=='/home/ruskin':
    poet_path=home_dir+'/projects/poet/'
    current_directory=home_dir+'/projects'+git_dir

sys.path.append(poet_path+'PythonPackage')
sys.path.append(poet_path+'scripts')

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

    parser.add_argument('-b',
                        action='store',
                        dest='breaks',
                        help='flag this to include tidal frequency breaks'
                        )

    return parser.parse_args()



if __name__ == '__main__':

    serialized_dir = poet_path +  "stellar_evolution_interpolators"
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    eccentricity_path=os.path.join(poet_path,'eccentricity_expansion_coef.txt').encode('ascii')

    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        eccentricity_path
    )


    args = cmdline_args()
    instance = args.instance
    system_number=args.system

    output_direcotry=current_directory+'/MCMC_'+system_number+'/'
    if os.path.isdir(output_direcotry)==False:os.mkdir(output_direcotry)

    catalog_file=current_directory+'/SpinlogQCatalog_el0.4.txt'
    stepfilename = output_direcotry+'step_file_'+args.instance+'.txt'

    with open(catalog_file,'r') as f:
        next(f)
        for lines in f:
            x=lines.split()
            at_system=x[0]
            if system_number==at_system:
                teff_value=float(x[2])
                teff_error=float(x[3])
                feh_value=float(x[4])
                feh_error=float(x[5])
                logg_value=float(x[10])
                logg_error=float(x[15])
                Porb_value=float(x[6])
                Porb_error=float(x[7])
                eccentricity_value=float(x[8])
                eccentricity_error=float(x[9])
                Pspin_value=float(x[12])
                Pspin_error=float(x[13])
                mass_ratio=float(x[14])
                break


    sampling_parameters=dict(age=dict(min=0.05,
                                      max=10.0,
                                      dist='Uniform',
                                      step=0.8),
                             primary_mass=dict(min=0.4,
                                               max=1.2,
                                               dist='Uniform',
                                               step=0.5),
                             feh=dict(value=feh_value,
                                      sigma=feh_error,
                                      dist='Normal',
                                      step=0.5),
                             Porb=dict(value=Porb_value,
                                       sigma=Porb_error,
                                       dist='Normal',
                                       step=Porb_error),
                             eccentricity=dict(value=eccentricity_value,
                                               sigma=eccentricity_error,
                                               dist='Normal',
                                               step=eccentricity_error),
                             Wdisk=dict(min=2*scipy.pi/14,
                                        max=2*scipy.pi/1.4,
                                        dist='Uniform',
                                        step=0.5),
                             logQ=dict(min=4.0,
                                       max=8.0,
                                       dist='Uniform',
                                       step=0.8)
                             )


    observational_parameters = dict(teff=dict(value=teff_value,
                                              sigma=teff_error),
                                    logg=dict(value=logg_value,
                                              sigma=logg_error),
                                    spin=dict(value=Pspin_value,
                                              sigma=Pspin_error
                                              )
                                    )

    fixed_parameters = dict(disk_dissipation_age=5e-3,
                            planet_formation_age=5e-3,
                            wind=True,
                            wind_saturation_frequency=2.54,
                            diff_rot_coupling_timescale=5e-3,
                            wind_strength=0.17,
                            inclination=scipy.pi/2
                            )

    print('Sampling Parameters: ',sampling_parameters)
    print('Observed Parameters: ',observational_parameters)


with open(stepfilename,'w') as f:
    for key,value in sampling_parameters.items():
        f.write(key + '\t' + repr(value['step']) + '\n')

mcmc = MetropolisHastings(system_number,
                          interpolator,
                          sampling_parameters,
                          fixed_parameters,
                          observational_parameters,
                          catalog_file,
                          mass_ratio,
                          instance,
                          output_direcotry)

sys.stdout.flush()
if args.start: mcmc.iterations()
elif args.cont: mcmc.continue_last()
else: print('provide correct arguments')

