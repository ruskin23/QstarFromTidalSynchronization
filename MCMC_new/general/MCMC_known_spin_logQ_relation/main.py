import argparse
import scipy
import sys
import os

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

    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        b"eccentricity_expansion_coef.txt"
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

    with open('data_file.txt','r') as f:
        for i,lines in enumerate(f):
            if i==data_line:
                data=lines.split()
                KIC=data[0]
                mass_ratio=float(data[13])
                break
    print(KIC)

    observation_data = dict(
                teff_primary=dict(value=float(data[1]),sigma=float(data[2])),
                feh=dict(value=float(data[3]),sigma=float(data[4])),
                Porb=dict(value=float(data[5]),sigma=float(data[6])),
                eccentricity=dict(value=float(data[7]),sigma=float(data[8])),
                logg=dict(value=float(data[9]),sigma=float(data[10])),
                        )

    observed_Pspin = dict(value=float(data[11]),sigma=float(data[12]))

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
                        teff_step=130.0,
                        feh_step=0.2,
                        #Porb_step=0.0001,
                        Porb_step=float(data[6]),
                        eccentricity_step=float(data[8]),
                        logg_step=0.1,
                        Wdisk_step=0.1,
                        logQ_step=0.15
                    )

    #proposed_step = dict(
                        #teff_step=10.0,
                        #feh_step=0.02,
                        #Porb_step=0.0001,
                        #eccentricity_step=float(data[6]),
                        #logg_step=0.01,
                        #Wdisk_step=0.1,
                        #logQ_step=0.15
                    #)

    Wdisk = dict(
                min=2*scipy.pi/14,
                max=2*scipy.pi/1.4
    )

    logQ = dict(value=float(data[14])
    )

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

