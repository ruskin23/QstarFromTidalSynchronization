import sys
from pathlib import Path
from directories import directories

home_dir=str(Path.home())
path=directories(home_dir)
sys.path.append(path.poet_path+'/PythonPackage')
sys.path.append(path.poet_path+'/scripts')
from orbital_evolution.transformations import phase_lag
from create_objects import BinaryObjects
import numpy


class IntialSecondaryAngmom:

    def __init__(self,
                 interpolator,
                 parameters
                 ):


        self.interpolator=interpolator
        self.parameters=parameters

        for item,value in parameters.items():
            setattr(self,item,value)

        if hasattr(parameters,'logQ'):self.convective_phase_lag=phase_lag(self.logQ)
        else: self.convective_phase_lag=self.phase_lag_max

    def __call__(self):

        binary_system=BinaryObjects(self.interpolator,self.parameters)

        star = binary_system.create_star(self.secondary_mass,dissipation=True)
        planet = binary_system.create_planet(1.0)

        binary = binary_system.create_binary_system(star,planet,initial_orbital_period=10.0,initial_eccentricity=0.0)

        binary.evolve(self.disk_dissipation_age, 1e-3, 1e-6, None)
        disk_state = binary.final_state()

        planet.delete()
        star.delete()
        binary.delete()

        return numpy.array([disk_state.envelope_angmom, disk_state.core_angmom])

