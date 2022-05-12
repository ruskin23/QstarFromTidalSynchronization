
import sys
import os
import random
import numpy
import scipy
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




serialized_dir =  "/home/ruskin/projects/poet/stellar_evolution_interpolators"
manager = StellarEvolutionManager(serialized_dir)
interpolator = manager.get_interpolator_by_name('default')

eccentricity_path=os.path.join('/home/ruskin/projects/poet','eccentricity_expansion_coef.txt').encode('ascii')

orbital_evolution_library.read_eccentricity_expansion_coefficients(
    eccentricity_path
)




interpolator('radius', 1.08939, -0.147).max_age