"""Allows finding mass given other age at fixed [Fe/H]."""
import sys
sys.path.append('/Users/ruskinpatel/Desktop/Research/poet/PythonPackage')
sys.path.append('/Users/ruskinpatel/Desktop/Research/poet/scripts')


from stellar_evolution.manager import StellarEvolutionManager
from stellar_evolution.derived_stellar_quantities import TeffK
import scipy
import scipy.interpolate
import scipy.linalg
import scipy.optimize

G=1534949.0910000005;

class DerivePrimnaryMass:

    def __init__(self, interpolator, feh, age, teff_value):
        """Set-up to use the given interpolator."""

        self.interp = interpolator
        self.feh = feh
        self.age = age
        self.teff_value = teff_value


    def teff_diff(self, mass):

        """Return the effective temperature for the given stellar mass."""


        T = TeffK(self.interp('radius', mass, self.feh),
                  self.interp('lum', mass, self.feh)
                 )


        if self.age >T.min_age and self.age<T.max_age:return (T(self.age)-self.teff_value)
        else : return scipy.nan

    def possible_solution(self):

        """Finds 2 values of masses between which a solution is present"""


        mass_array = scipy.linspace(0.4,1.2,100)
        teff_array_diff = scipy.empty(len(mass_array))

        for m_index,m_value in enumerate(mass_array):

            teff_array_diff[m_index] = self.teff_diff(m_value)

        mass_solutions = []

        for i in range(teff_array_diff.size-1):

            x = teff_array_diff[i]*teff_array_diff[i+1]
            if x<0:
                mass_solutions.append(mass_array[i])
                mass_solutions.append(mass_array[i+1])
                return mass_solutions
                break


        if i >= teff_array_diff.size-1 : return scipy.nan

    def __call__(self):

        """Optimize the possible solution to calculate exact solution"""

        mass_solutions = self.possible_solution()

        solution = scipy.optimize.brentq(self.teff_diff, mass_solutions[0], mass_solutions[1])

        return solution


serialized_dir = "/home/kpenev/projects/git/poet/stellar_evolution_interpolators"
manager = StellarEvolutionManager(serialized_dir)
interpolator = manager.get_interpolator_by_name('default')

x=DerivePrimnaryMass(interpolator,-0.0329537979168565, 12.918385368071597 ,  5643.472584167178)

print(interpolator('radius', 1.0, -0.1))
#print(x())
