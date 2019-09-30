import argparse
import scipy
import sys
import os
import os.path

from pathlib import Path
home_dir=str(Path.home())

if home_dir=='/home/rxp163130':
    poet_path=home_dir+'/poet/'
if home_dir=='/home/ruskin':
    poet_path=home_dir+'/projects/poet/'

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

    args = parser.parse_args()

    system_number=args.system
    data_line=int(system_number)
    data_filename = os.getcwd() + '/data_file.txt'
    with open(data_filename,'r') as f:
        for i,lines in enumerate(f):
            if i==data_line:
                data=lines.split()
                KIC=data[0]
                mass_ratio=float(data[14])
                break
    print(KIC)
    print(mass_ratio)
    mass_age_feh_sample_file=os.getcwd() + '/mass_age_feh_sample_'+system_number+'.txt'

    observation_data = dict(
                Porb=dict(value=float(data[6]),sigma=float(data[7])),
                eccentricity=dict(value=float(data[8]),sigma=float(data[9])),
                        )

    print(observation_data)

    observed_Pspin = dict(value=float(data[12]),sigma=float(data[13]))

    print(observed_Pspin)

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
                        Porb_step=float(data[7]),
                        eccentricity_step=float(data[9]),
                        Wdisk_step=0.1,
                        logQ_step=1.0
                    )


    Wdisk = dict(
                min=2*scipy.pi/14,
                max=2*scipy.pi/1.4
    )

    logQ = dict(value=float(data[15])
    )

    print(logQ)
output_direcotry= os.getcwd()+'/MCMC_'+system_number+'/'
if os.path.isdir(output_direcotry)==False:os.mkdir(output_direcotry)

stepfilename = output_direcotry+'step_file_'+args.instance+'.txt'
with open(stepfilename,'w') as f:
    for key,value in proposed_step.items():
        f.write(key + '\t' + repr(value) + '\n')

instance = args.instance

mcmc = MetropolisHastings(
                            interpolator,
                            fixed_parameters,
                            observation_data,
                            mass_age_feh_sample_file,
                            Wdisk,
                            logQ,
                            proposed_step,
                            observed_Pspin,
                            mass_ratio,
                            system_number,
                            instance,
                            output_direcotry)


if args.start: mcmc.iterations()
elif args.cont: mcmc.continue_last()
else: print('provide correct arguments')

