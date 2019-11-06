import argparse
import scipy
import sys
import os
import os.path

from pathlib import Path
home_dir=str(Path.home())

git_dir='/QstarFromTidalSynchronization/MCMC/general/'
samples_dir='/QstarFromTidalSynchronization/MCMC/mcmc_mass_age/samples'

if home_dir=='/home/rxp163130':
    poet_path=home_dir+'/poet/'
    current_directory=home_dir+git_dir
    samples_directory=home_dir+samples_dir

if home_dir=='/home/ruskin':
    poet_path=home_dir+'/projects/poet/'
    current_directory=home_dir+'/projects'+git_dir
    samples_directory=home_dir+'/projects'+samples_dir


sys.path.append(poet_path+'PythonPackage')
sys.path.append(poet_path+'scripts')

from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from metropolis_hasting import MetropolisHastings


if __name__ == '__main__':

    serialized_dir = poet_path +  "stellar_evolution_interpolators"
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    eccentricity_path=os.path.join(poet_path,'eccentricity_expansion_coef.txt').encode('ascii')

    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        eccentricity_path
    )

    parser = argparse.ArgumentParser()
    parser.add_argument('-s', action='store_const', dest='start',
                    const='start',
                    help='start mcmc from beginning')

    parser.add_argument('-c', action='store_const', dest='cont',
                    const='continue',
                    help='continue mcmc from last iteration')

    parser.add_argument('-i', action = 'store', dest = 'instance',
                    help = 'define an instance of mcmc')

    parser.add_argument('-l', action = 'store', dest = 'system',
                    help = 'select a system for mcmc')

    parser.add_argument('-b', action = 'store', dest = 'breaks',
                    help = 'flag this to include tidal frequency breaks')

    args = parser.parse_args()

    system_number=args.system

    mass_age_feh_sample_file=samples_directory+ '/MassAgeFehSamples_'+system_number+'.txt'
    catalog_file=current_directory+'/spin_vs_logQ_systems_0.2.txt'
    solution_file=current_directory+'/SolutionFile.txt'


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

    observation_data = dict(Porb=dict(value=Porb_value,sigma=Porb_error),
                            eccentricity=dict(value=eccentricity_value,sigma=eccentricity_error)
                        )

    print('Observational data: ',observation_data)

    observed_Pspin = dict(value=Pspin_value,sigma=Pspin_error)

    print('Observed Spin: ',observed_Pspin)

    fixed_parameters = dict(
                        disk_dissipation_age=5e-3,
                        planet_formation_age=5e-3,
                        wind=True,
                        wind_saturation_frequency=2.54,
                        diff_rot_coupling_timescale=5e-3,
                        wind_strength=0.17,
                        inclination=scipy.pi/2

    )

    proposed_step = dict(
                        Porb_step=Porb_error,
                        eccentricity_step=eccentricity_error,
                        Wdisk_step=0.1,
                        logQ_step=0.2
                    )


    Wdisk = dict(
                min=2*scipy.pi/14,
                max=2*scipy.pi/1.4
    )

    with open(solution_file,'r') as f:
        next(f)
        for lines in f:
            x=lines.split()
            at_system=x[0]
            if at_system==system_number:
                logQ_value=float(x[1])
                break
    logQ = dict(value=logQ_value)

    print('logQ: ',logQ)

output_direcotry=current_directory+'/MCMC_'+system_number+'/'

if os.path.isdir(output_direcotry)==False:os.mkdir(output_direcotry)

stepfilename = output_direcotry+'step_file_'+args.instance+'.txt'
with open(stepfilename,'w') as f:
    for key,value in proposed_step.items():
        f.write(key + '\t' + repr(value) + '\n')

instance = args.instance

mcmc = MetropolisHastings(system_number,
                          interpolator,
                          fixed_parameters,
                          observation_data,
                          mass_age_feh_sample_file,
                          catalog_file,
                          Wdisk,
                          logQ,
                          proposed_step,
                          observed_Pspin,
                          mass_ratio,
                          instance,
                          output_direcotry)


if args.start: mcmc.iterations()
elif args.cont: mcmc.continue_last()
else: print('provide correct arguments')

