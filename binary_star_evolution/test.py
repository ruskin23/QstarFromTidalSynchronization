import sys
sys.path.append('.../poet/PythonPackage')
sys.path.append('.../poet/scripts')

from matplotlib import pyplot
from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as \
    orbital_evolution_library
from orbital_evolution.binary import Binary
from orbital_evolution.transformations import phase_lag
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet
from basic_utils import Structure
from astropy import units, constants
from basic_utils import Structure
from math import pi
from scipy.optimize import brentq
import numpy







if __name__ == '__main__':
    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        b"eccentricity_expansion_coef.txt"
    )
    serialized_dir = '/home/kpenev/projects/git/poet/stellar_evolution_interpolators'
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    # test_evolution(interpolator, phase_lag(6.0))
    logQ = numpy.linspace(5.5,6.0,10)

    with open ('test.txt', 'a') as f:


        for q in logQ:
        
            f.write("Q = " + "\t" +  repr(phase_lag(q)) + "\tlogvalue" + repr(q) + "\n")


