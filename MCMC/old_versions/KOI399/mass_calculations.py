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
import scipy.interpolate
import scipy.linalg
import scipy.optimize

G = 1534949.0910000005;
constant = 0.0244e7


class DeriveMass:

    def __init__(self, interpolator, feh, age, teff_value):
        """Set-up to use the given interpolator."""

        self.interp = interpolator
        self.feh = feh
        self.age = age
        self.teff_value = teff_value

        self.mass_bound_check = False

    def teff_diff(self, mass):

        """Return the effective temperature for the given stellar mass."""



        T = TeffK(self.interp('radius', mass, self.feh),
                  self.interp('lum', mass, self.feh)
                  )

        if self.age >T.min_age and self.age<T.max_age:

            return (T(self.age) - self.teff_value)

        else: return scipy.nan

    def possible_solution(self):

        """Finds 2 values of masses between which a solution is present"""

        mass_array = scipy.linspace(0.4, 1.2, 100)
        teff_array_diff = scipy.empty(len(mass_array))

        for m_index, m_value in enumerate(mass_array):

            teff_array_diff[m_index] = self.teff_diff(m_value)


        mass_solutions = []
        os.remove('mass_test.txt')
        for i in range(teff_array_diff.size - 1):
            x = teff_array_diff[i] * teff_array_diff[i + 1]
            with open('mass_test.txt', 'a') as f:
                f.write(repr(i) + '\t' + repr(x) + '\n')
            if x < 0:
                mass_solutions.append(mass_array[i])
                mass_solutions.append(mass_array[i + 1])
                return mass_solutions

        if i >=teff_array_diff.size - 2: self.mass_bound_check = True

    def __call__(self):

        """Optimize the possible solution to calculate exact solution"""

        solution = 0
        mass_solutions = self.possible_solution()
        print(mass_solutions)

        if (self.mass_bound_check) is False:
            solution = scipy.optimize.brentq(self.teff_diff, mass_solutions[0], mass_solutions[1])
            print("solution = ", solution)
            return solution
        else: return scipy.nan


