#!/usr/bin/env python3 -u

import argparse
import scipy
import sys
import os
import os.path

from pathlib import Path
home_dir=str(Path.home())

poet_path=home_dir+'/projects/poet/'
current_directory=home_dir+'/projects'+git_dir
samples_directory=home_dir+'/projects/QstarFromTidalSynchronization/MCMC/mcmc_mass_age/samples/updated_samples'

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

    parser.add_argument('-m',
                        action='store',
                        dest='sampling_method',
                        help='select sampling method:uncorrelated or adaptive')


    parser.add_argument('-t',
                        action-'store',
                        dest='test_case',
                        help='choose from following cases:
                        prior:L=1.0
                        gp:L=N(phi;sigma_phi)
                        gpt:L=N(phi;sigma_phi)N(theta;sigma_theta)
                        gptc:L=N(phi,theta;sigma_phi_theta)')

    parser.add_argument('-b',
                        action='store',
                        dest='breaks',
                        help='flag this to include tidal frequency breaks'
                        )

    return parser.parse_args()



if __name__ == '__main__':

    args = cmdline_args()
    instance = args.instance
    system_number=args.system

    test_case=args.test_case
    sampling_method=args.sampling_method
    if sampling_method=='uncorrelated':
        sys.path.append('../Uncorrelated_Phi')
        output_directory='../Uncorrelated_Phi/'+test_case
    if sampling_method=='adaptive':
        sys.path.append('../Adaptive')
        output_directory='../Adaptive/'+test_case
    if os.path.isdir(output_directory)==False:os.mkdir(output_directory)

    catalog_file=current_directory+'/SpinlogQCatalog_el0.4.txt'
    solution_file=current_directory+'/SolutionFileBreaks0.0.txt'
    samples_file=samples_directory+'/MassAgeFehSamples_'+system_number+'.txt'

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


    print('Sampling Parameters: ',sampling_parameters)


mcmc = MetropolisHastings(system_number,
                          interpolator,
                          sampling_parameters,
                          mass_ratio,
                          instance,
                          sampling_method,
                          test_case,
                          output_directory)

sys.stdout.flush()
if args.start: mcmc.iterations()
elif args.cont: mcmc.continue_last()
else: print('provide correct arguments')

