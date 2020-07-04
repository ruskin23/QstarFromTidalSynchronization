#!/usr/bin/env python3 -u

import argparse
import scipy

import os
import sys

poet_path='/home/ruskin/projects/poet/'
sys.path.append(poet_path+'PythonPackage')
sys.path.append(poet_path+'scripts')

from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from sampling_test_class import NestedSampling


def cmdline_args():

    parser = argparse.ArgumentParser()

    parser.add_argument('-l',
                        action='store',
                        dest='system',
                        help='select a system for mcmc'
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
    system_number=args.system

    catalog_file='SpinlogQCatalog_el0.4.txt'

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


    sampling_parameters=dict(feh=dict(value=feh_value,
                                      sigma=feh_error,
                                      low=-1.014,
                                      high=0.537,
                                      dist='Truncated_Normal'),
                             Porb=dict(value=Porb_value,
                                       sigma=Porb_error,
                                       dist='Normal'),
                             eccentricity=dict(value=eccentricity_value,
                                               sigma=eccentricity_error,
                                               low=0.0,
                                               high=0.45,
                                               dist='Truncated_Normal'),
                             Wdisk=dict(low=2*scipy.pi/14,
                                        high=2*scipy.pi/1.4,
                                        dist='Uniform'),
                             logQ=dict(low=5.0,
                                       high=12.0,
                                       dist='Uniform'),
                             primary_mass=dict(low=0.5,
                                               high=1.2,
                                               dist='Uniform'),
                             age=dict(low=1e-3,
                                      high=10,
                                      dist='Uniform')

                             )


    observed_parameters=dict(teff=dict(value=teff_value,
                                       sigma=teff_error),
                             logg=dict(value=logg_value,
                                       sigma=logg_error),
                             spin=dict(value=Pspin_value,
                                       sigma=Pspin_error)
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
    print('Observed Parameters: ',observed_parameters)

sampling = NestedSampling(system_number,
                          interpolator,
                          sampling_parameters,
                          fixed_parameters,
                          observed_parameters,
                          mass_ratio)


sampling.start()
