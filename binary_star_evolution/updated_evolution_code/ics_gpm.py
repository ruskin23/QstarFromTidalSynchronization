#!/usr/bin/env python3

from copyreg import pickle
import logging
import sys
from pathlib import Path
from directories import directories

home_dir=str(Path.home())
path=directories(home_dir)
sys.path.append(path.poet_path+'/PythonPackage')
sys.path.append(path.poet_path+'/scripts')
sys.path.append('/home/ruskin/projects')


from orbital_evolution.transformations import phase_lag
from create_objects import BinaryObjects

import numpy
from astropy import units
from types import SimpleNamespace

from general_purpose_python_modules.reproduce_system import find_evolution

_logger = logging.getLogger(__name__)


class InitialConditionSolver:
    """Find initial orbital period and eccentricity which reproduce
        current orbital period and eccentricity of a given system """


    def __init__(self,
                interpolator,
                parameters,
                evolution_max_time_step=1e-3,
                evolution_precision=1e-6,
                secondary_angmom=None,
                initial_guess=None):
        """
        Args:

            - interpolator

            - parameters
                a dictionary of all the parameters

            - evolution_max_time_step:
                The maximum timestep the evolution is allowed to make.

            - evolution_precision:
                The precision to require of the evolution.

            - secondary_angmom:
                The initial angular momentum of secondary star
        """

        self.interpolator=interpolator
        self.parameters=parameters
        for item,value in parameters.items():
            setattr(self,item,value)
        if 'logQ' in parameters: self.convective_phase_lag=phase_lag(self.logQ)
        else: self.convective_phase_lag=self.phase_lag_max

        self.evolution_max_time_step=evolution_max_time_step
        self.evolution_precision = evolution_precision
        self.secondary_angmom = secondary_angmom


    def get_dissipation_dict(self, lgq_suffixes=('primary', 'secondary')):

        result = dict()
        for component in lgq_suffixes:
            
            result[component] = dict()
            result[component]['reference_phase_lag'] = self.convective_phase_lag
            result[component]['spin_frequency_breaks']=None
            result[component]['spin_frequency_powers']=numpy.array([0.0])
            result[component]['tidal_frequency_breaks']=numpy.array([0.12566371,0.30692368])
            result[component]['tidal_frequency_powers']=numpy.array([ 0.0,0.19434417,0.0])

        return result


    def __call__(self, primary, secondary):
        """
        Find initial conditions which reproduce the given system now.

        Args:
            - primary

            - secondary
        """

        self.primary = primary
        self.secondary = secondary

        self.target_age=self.age
        self.target_orbital_period=self.orbital_period
        self.target_eccentricity=self.eccentricity
        
        
        kwargs = dict(system=SimpleNamespace(primary_mass=self.primary_mass * units.M_sun,
                                             secondary_mass=self.secondary_mass * units.M_sun,
                                             feh=self.feh,
                                             orbital_period=self.target_orbital_period* units.day,
                                             eccentricity=self.target_eccentricity,
                                             age=self.age* units.Gyr),
                    interpolator=self.interpolator,
                    dissipation=self.get_dissipation_dict(),
                    max_age=self.age* units.Gyr,
                    disk_period=2*numpy.pi/self.Wdisk* units.day,
                    disk_dissipation_age=self.disk_dissipation_age* units.Gyr,
                    primary_core_envelope_coupling_timescale=self.diff_rot_coupling_timescale* units.Gyr,
                    secondary_core_envelope_coupling_timescale=self.diff_rot_coupling_timescale* units.Gyr,
                    secondary_wind_strength=self.wind_strength,
                    primary_wind_saturation=self.wind_strength,
                    secondary_wind_saturation=self.wind_saturation_frequency,
                    required_ages=None,
                    timeout=3600,
                    initial_eccentricity=True,
                    orbital_period_tolerance=1e-4,
                    period_search_factor=2,
                    scaled_period_guess=2,
                    max_time_step=self.evolution_max_time_step,
                    precision=self.evolution_precision)

        evolution=find_evolution(**kwargs)
        with open('reproduced_evolution.pickle','rb') as f:
            pickle.dump(evolution,f)