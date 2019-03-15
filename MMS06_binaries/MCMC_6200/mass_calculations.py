"""Calculate primary and secondary mass if the age, [Fe/H] and semimajor axis are known """


import sys
import os

#sys.path.append('/Users/ruskinpatel/Desktop/Research/poet/PythonPackage')
#sys.path.append('/Users/ruskinpatel/Desktop/Research/poet/scripts')

sys.path.append('/home/kpenev/projects/git/poet/PythonPackage')
sys.path.append('/home/kpenev/projects/git/poet/scripts')

from stellar_evolution.manager import StellarEvolutionManager
from stellar_evolution.derived_stellar_quantities import TeffK
import scipy
import numpy
import scipy.optimize
import matplotlib.pyplot as plt

G = 15.3588801038e5

class DeriveSecondaryMass:

    def __init__(self, porb, vkr, i, mass_primary):

        self.porb = porb
        self.vkr = vkr
        self.i = i
        self.mass_primary = mass_primary
        self.mass_found = None

    def difference_function(self, mass_secondary):

        mass_function = (mass_secondary*scipy.sin(self.i))**3 /(mass_secondary + self.mass_primary)**2

        return mass_function - (self.porb * (self.vkr)**3)/(2*scipy.pi*G)

    def __call__(self):

        mass_array = scipy.linspace(0.005,1.4,100)
        differnece_array = scipy.empty(len(mass_array))

        for m_index, m_value in enumerate(mass_array):
            differnece_array[m_index] = self.difference_function(m_value)

        mass_solutions = []

        for i in range(differnece_array.size-1):

            x = differnece_array[i]*differnece_array[i+1]
            if x<0:
                mass_solutions.append(mass_array[i])
                mass_solutions.append(mass_array[i+1])
                self.mass_found=True
                break

        if self.mass_found:
            solution = scipy.optimize.brentq(self.difference_function, mass_solutions[0], mass_solutions[1])

            print('PRIMARY MASS = ', self.mass_primary)
            print ('SECONDARY MASS = ', solution)

            return solution

        else : return scipy.nan



