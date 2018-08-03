"""Allows finding mass given other age at fixed [Fe/H]."""
import sys
sys.path.append('..')


from stellar_evolution.library_interface import MESAInterpolator
from stellar_evolution.derived_stellar_quantities import\
    TeffK,\
    LogGCGS,\
    RhoCGS
from basic_utils import Structure
import scipy
import scipy.interpolate
import scipy.linalg
import scipy.optimize

class DerivePrimnaryMass:

    def __init__(self, interpolator, feh, age, teff_value):
        """Set-up to use the given interpolator."""

        self.interp = interpolator
        self.feh = feh
        self.age = age
        self.teff_value = teff_value


    def teff(self, mass):

        """Return the effective temperature for the given stellar mass."""

        T = TeffK(self.interp('radius', mass, self.feh),
                  self.interp('lum', mass, self.feh)
                 )

        return T(self.age)


    def possible_solution(self,value):

        """Finds 2 values of masses between which a solution is present"""


        mass_array = scipy.linspace(0.01,2.0,100)
        teff_array = scipy.empty(len(mass_array))
        teff_array_diff = scipy.empty(len(mass_array))

        for m_index,m_value in enumerate(mass_array):

            teff_array[m_index] = self.teff(m_value)

        teff_array_diff[:] = [x - value for x in t_array]

        mass_solutions = scipy.empty(2)

        for i in teff_array_diff:

            x = teff_array_diff[i]*teff_array_diff[i+1]
            if x<0:
                mass_solutions.append(mass_array[i])
                mass_solutions.append(mass_array[i+1])
                break


        return mass_solutions


    def find_solution(self,value):

        """Optimize the possible solution to calculate exact solution"""

        mass_solutions = self.possible_solution(value)

        solution = scipy.optimize.brentq(self.teff,mass_solutions[0],mass_solutions[1])

        return solution













