"""Allows finding mass given other age at fixed [Fe/H]."""
import sys
sys.path.append('..')


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

        return (T(self.age)-self.teff_value)


    def possible_solution(self):

        """Finds 2 values of masses between which a solution is present"""


        mass_array = scipy.linspace(0.4,1.2,100)
        teff_array_diff = scipy.empty(len(mass_array))

        for m_index,m_value in enumerate(mass_array):

            teff_array_diff[m_index] = self.teff_diff(m_value)

        mass_solutions = []

        for i in range(teff_array_diff.size):

            x = teff_array_diff[i]*teff_array_diff[i+1]
            if x<0:
                mass_solutions.append(mass_array[i])
                mass_solutions.append(mass_array[i+1])
                break


        return mass_solutions


    def __call__(self):

        """Optimize the possible solution to calculate exact solution"""

        mass_solutions = self.possible_solution()

        solution = scipy.optimize.brentq(self.teff_diff, mass_solutions[0], mass_solutions[1])

        return solution


class DeriveSecondaryMass:


    def __init__(self, p_orb, v_r, i):
	
        self.p_orb = p_orb
        self.v_r = v_r
        self.i = i

    def mass_function(self,mass_secondary,mass_primary):

        return (mass_secondary**3 * scipy.sin(self.i)**3 )/(mass_secondary + self.mass_primary)**2

    def difference_function(self, mass_secondary,mass_primary):

        return self.mass_function(mass_secondary,mass_primary) - (self.p_orb * (self.v_r * scipy.sin(self.i))**3)/2*scipy.pi*G

    def __call__(self,mass_primary):

        mass_array = scipy.linspace(0.4,1.2,100)
        differnece_array = scipy.empty(len(mass_array))

        for m_index, m_value in enumerate(mass_array):
            differnece_array[m_index] = self.difference_function(m_value,mass_primary)

        mass_solutions = []

        for i in range(differnece_array.size):

            x = differnece_array[i]*differnece_array[i+1]
            if x<0:
                mass_solutions.append(mass_array[i])
                mass_solutions.append(mass_array[i+1])
                break

        solution = scipy.optimize.brentq(self.difference_function, mass_solutions[0], mass_solutions[1])

        return solution

serialized_dir = "/home/kpenev/projects/git/poet/stellar_evolution_interpolators"
manager = StellarEvolutionManager(serialized_dir)
interpolator = manager.get_interpolator_by_name('default')

x=DerivePrimnaryMass(interpolator, 0, 4.56, 5870)
primary_mass = x.find_solution()
