import argparse
import os
import sys
from pathlib import Path
from directories import directories

home_dir=str(Path.home())
path=directories(home_dir)
sys.path.append(path.poet_path+'/PythonPackage')
sys.path.append(path.poet_path+'/scripts')

from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from stellar_evolution.derived_stellar_quantities import\
    TeffK,\
    LogGCGS,\
    RhoCGS

import scipy
from scipy import stats
import numpy
import logging
from test_spin_calculation import SpinPeriod
import time
from datetime import date
from pathos.pools import ProcessPool

today=date.today()
today=today.strftime("%B_%d")

def cmdline_args():

    parser = argparse.ArgumentParser()

    parser.add_argument('-s',
                        action='store',
                        dest='status',
                        help='starting or continuing run')

    parser.add_argument('-l',
                        action='store',
                        dest='system',
                        help='select a system for mcmc'
                        )

    parser.add_argument('-n',
                        action='store',
                        dest='nprocs',
                        help='number threads')


    parser.add_argument('-i',
                        action='store',
                        dest='instance',
                        help='instance')

    return parser.parse_args()

def likelihood_logger(scratch_directory,system):
    l_logger=logging.getLogger('spin_calculations')
    if not l_logger.handlers:
        l_logger.setLevel(logging.DEBUG)
        pid=str(os.getpid())[-3:]
        scratch_filename=f'{scratch_directory}/system_{system}/test_{today}_sampling_{pid}.log'
        proces_handler=logging.FileHandler(scratch_filename)

        process_format=logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s',datefmt='%d-%b-%y %H:%M:%S')
        proces_handler.setFormatter(process_format)

        l_logger.addHandler(proces_handler)
    return l_logger


def calculate_model(PARAMS):

    system=PARAMS[0]
    interpolator=PARAMS[1]
    parameter_set=PARAMS[2]
    fixed_parameters=PARAMS[3]
    scratch_directory=PARAMS[4]
    sampling_logger=likelihood_logger(scratch_directory,system)

    model_set=dict()

    sampling_logger.info('Parameter Set')
    print(parameter_set)
    for key,value in parameter_set.items():
        sampling_logger.info('{} = {}'.format(key,value))
        print('{} = {}'.format(key,value))


    mass=parameter_set['primary_mass']
    feh=parameter_set['feh']
    age=parameter_set['age']

    quantity_radius=interpolator('radius',mass,feh)
    quantity_lum=interpolator('lum',mass,feh)

    T=TeffK(quantity_radius,quantity_lum)
    try:model_set['teff']=T(age)
    except:model_set['teff']=-scipy.inf
    G=LogGCGS(mass,quantity_radius)
    try:model_set['logg']=G(age)
    except:model_set['logg']=-scipy.inf

    spin_calculations = SpinPeriod(
                                system,
                                interpolator,
                                parameter_set,
                                fixed_parameters,
                                mass_ratio,
                                sampling_logger
                                )

    model_set['spin']=spin_calculations()
    if numpy.isnan(model_set['spin']):model_set['spin']=-scipy.inf
    return model_set['spin']

        # return model_set


if __name__ == '__main__':

    args = cmdline_args()
    status=args.status
    system_number=args.system
    instance=args.instance
    nprocs=int(args.nprocs)
    if instance is not None:instance=int(instance)


    serialized_dir = path.poet_path +  "/stellar_evolution_interpolators"
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    eccentricity_path=os.path.join(path.poet_path,'eccentricity_expansion_coef.txt').encode('ascii')

    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        eccentricity_path
    )



    catalog_file=path.current_directory+'/SpinlogQCatalog_el0.4.txt'

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


    sampling_parameters = [('Porb',10.303032034625893,Porb_error,'Normal'),
                           ('feh', -0.31882201216054464,feh_error,-1.014,0.537,'Turncated_Normal'),
                           ('eccentricity', 0.27557873381900033,eccentricity_error,0.0,0.45,'Turncated_Normal'),
                           ('Wdisk',  1.786126160343231,2*scipy.pi/14,2*scipy.pi/1.4,'Uniform'),
                           ('logQ',  5.6832393389266604,5.0,12.0,'Uniform'),
                           ('primary_mass',  0.7123097489857166,0.5,1.2,'Uniform'),
                           ('age',  4.525614334426742,1e-3,10.0,'Uniform')]


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
    print(f'Sampling Parameters: {sampling_parameters}')
    print(f'Observed Parameters: {observed_parameters}')

    scratch_directory=path.scratch_directory

    parameter_set=dict()
    parameter_set_1=dict()
    parameter_set_2=dict()
    parameter_set_3=dict()
    parameter_set_4=dict()


    for s in sampling_parameters:
        parameter_set[s[0]]=s[1]
        parameter_set_1[s[0]]=scipy.stats.norm.rvs(loc=s[1],scale=0.01)
        parameter_set_2[s[0]]=scipy.stats.norm.rvs(loc=s[1],scale=0.01)
        parameter_set_3[s[0]]=scipy.stats.norm.rvs(loc=s[1],scale=0.01)
        parameter_set_4[s[0]]=scipy.stats.norm.rvs(loc=s[1],scale=0.01)
    

    param=(system_number,interpolator,parameter_set,fixed_parameters,scratch_directory)
    param_1=(system_number,interpolator,parameter_set_1,fixed_parameters,scratch_directory)
    param_2=(system_number,interpolator,parameter_set_2,fixed_parameters,scratch_directory)
    param_3=(system_number,interpolator,parameter_set_3,fixed_parameters,scratch_directory)
    param_4=(system_number,interpolator,parameter_set_4,fixed_parameters,scratch_directory)

    PARAMS=[param,param_1,param_2,param_3,param_4]

    pool=ProcessPool(nodes=nprocs)

    # p=numpy.array(list(pool.map(calculate_model,PARAMS)))
    pool.map(calculate_model,PARAMS)
    

    # for i in range(3):
    
    #     sampling_logger=likelihood_logger(scratch_directory,system_number)
    #     model_set=calculate_model(system_number,interpolator,parameter_set,fixed_parameters,sampling_logger)
    #     L=1
    #     sampling_logger.info(f'Loglike = {L}')

    #     for key,value in parameter_set.items():
    #         parameter_set[key]=scipy.stats.norm.rvs(loc=parameter_set[key],scale=0.01)
        

        
