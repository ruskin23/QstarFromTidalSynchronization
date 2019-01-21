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
    logQ = numpy.linspace(5.5,6.0,10)
    x = numpy.linspace(0,12,1000)
    p_value = 7.713253717543052
    with open ('current_psin.txt', 'a') as f:


        for q in x:
        
            f.write(repr(q) + "\t" + repr(p_value) + "\n")


