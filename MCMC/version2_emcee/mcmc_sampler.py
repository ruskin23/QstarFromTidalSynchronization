import sys
import os
import random
import numpy
sys.path.append('/home/ruskin/projects/poet/PythonPackage')
sys.path.append('/home/ruskin/projects/poet/scripts')
sys.path.append('/home/ruskin/projects/QstarFromTidalSynchronization/binary_star_evolution/updated_evolution_code')

from stellar_evolution.library_interface import library
from create_objects import BinaryObjects
from prior_transform_class import prior_transform
from initial_condition_solver import InitialConditionSolver
from initial_secondary_angmom import IntialSecondaryAngmom

from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as \
    orbital_evolution_library
from orbital_evolution.transformations import phase_lag


class sampler:

    def __init__(self,
                system_number,
                uniform_variables):

        self.system_number=system_number
        self.uniform_variable=uniform_variables

        self.params=dict()

    def param_conversion(self,
                        sampled_params):

        sum_mass=sampled_params[0]
        mass_ratio=sampled_params[1]

        primary_mass=(1/(1+mass_ratio))*sum_mass
        secondary_mass=primary_mass*mass_ratio

        return [primary_mass,secondary_mass]

    def sampled_from_data(self):

        params=prior_transform(self.system_number)
        sampled_params=params.paramter_evaluate(self.uniform_variable)
        masses=self.param_conversion(sampled_params)


        self.params['primary_mass'] = masses[0]
        self.params['secondary_mass'] = masses[1]
        z=sampled_params[2]
        self.params['feh'] = library.feh_from_z(z)
        self.params['age'] = sampled_params[3]
        self.params['eccentricity'] = sampled_params[4]

    def fixed_params(self):

        self.params['dissipation']=True
        self.params['disk_dissipation_age']=5e-3
        self.params['wind']=True
        self.params['wind_saturation_frequency']=2.54
        self.params['diff_rot_coupling_timescale']=5e-3
        self.params['wind_strength']=0.17

    def uniformly_sampled(self):

        phase_lag_max=numpy.random.uniform(phase_lag(5),phase_lag(12))
        # phase_lag_max=phase_lag(5)
        alpha = numpy.random.uniform(1,5)

        lnP= numpy.random.uniform(numpy.log(0.5),numpy.log(50))
        omegaref= numpy.exp(numpy.log(2*numpy.pi) - lnP) 
        omegamin=(2*numpy.pi)/50

        if alpha<0:tidal_power=numpy.array([1,0,alpha])
        else:tidal_power=numpy.array([1,alpha,0])

        self.params['phase_lag_max']=phase_lag_max

        self.params['tidal_frequency_breaks']=numpy.array([omegamin,omegaref])
        self.params['tidal_frequency_powers']=tidal_power


        # self.params['tidal_frequency_breaks']=None
        # self.params['tidal_frequency_powers']=numpy.array([0.0])
        self.params['spin_frequency_breaks']=None
        self.params['spin_frequency_powers']=numpy.array([0.0])


        self.params['Wdisk'] = numpy.random.uniform(2*numpy.pi/14,2*numpy.pi/1.4)

    def get_orbital_period(self):
        with open('/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/version2_emcee/catalog/filtering/orbital_periods.txt','r') as f:
            for lines in f:
                x=lines.split()
                if self.system_number==x[0]:
                    self.params['orbital_period']=float(x[1])
                    break


    def __call__(self):

        self.uniformly_sampled()
        self.fixed_params()
        self.sampled_from_data()
        self.get_orbital_period()

        return self.params





if __name__ == '__main__':

    serialized_dir =  "/home/ruskin/projects/poet/stellar_evolution_interpolators"
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    eccentricity_path=os.path.join('/home/ruskin/projects/poet','eccentricity_expansion_coef.txt').encode('ascii')

    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        eccentricity_path
    )



    system_number='10296163'
    
    while(True):
        uniform_variables=[random.uniform(0,1),random.uniform(0,1),random.uniform(0,1),random.uniform(0,1),random.uniform(0,1)]
        get_params=sampler(system_number,uniform_variables)
        parameters=get_params()

        if numpy.logical_or(parameters['feh']<-1.014,
                            parameters['feh']>0.537):
                            continue
        
        if numpy.logical_or(numpy.logical_or(parameters['primary_mass']<0.4,
                                             parameters['primary_mass']>1.2),
                            numpy.logical_or(parameters['secondary_mass']<0.4,
                                            parameters['secondary_mass']>1.2)
                            ):
                            continue

        for key,values in parameters.items():
            print(f'{key}\t{values}\n')

        stars=BinaryObjects(interpolator,parameters)
        primary=stars.create_star(parameters['primary_mass'])
        secondary=stars.create_star(parameters['secondary_mass'])
        
        angmom=IntialSecondaryAngmom(interpolator,parameters)


        initial_conditions=InitialConditionSolver(interpolator,parameters,secondary_angmom=angmom())
        try:
            results=initial_conditions(primary,secondary)
        except:
            continue

        if numpy.isnan(results['spin']):
            print('\n')
            continue
        else: break