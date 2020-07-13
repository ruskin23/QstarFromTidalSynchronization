"""Calculate primary and secondary mass if the age, [Fe/H] and semimajor axis are known """


import sys
import os

import os.path

from pathlib import Path
home_dir=str(Path.home())

git_dir='/QstarFromTidalSynchronization/MCMC/combined'

if home_dir=='/home/rxp163130':
    poet_path=home_dir+'/poet/'
    current_directory=home_dir+git_dir
    samples_directory=home_dir+'/QstarFromTidalSynchronization/MCMC/mcmc_mass_age/samples/updated_samples'

if home_dir=='/home/ruskin':
    poet_path=home_dir+'/projects/poet/'
    current_directory=home_dir+'/projects'+git_dir
    samples_directory=home_dir+'/projects/QstarFromTidalSynchronization/MCMC/mcmc_mass_age/samples/updated_samples'

if home_dir=='/home1/06850/rpatel23':
    work_dir='/work/06850/rpatel23/stampede2'
    poet_path=work_dir+'/poet'
    current_directory=work_dir+git_dir
    samples_directory=work_dir+'/QstarFromTidalSynchronization/MCMC/mcmc_mass_age/samples/updated_samples'

sys.path.append(poet_path+'PythonPackage')
sys.path.append(poet_path+'scripts')

from stellar_evolution.manager import StellarEvolutionManager
from stellar_evolution.derived_stellar_quantities import TeffK
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
import scipy
import numpy
import scipy.interpolate
import scipy.linalg
import scipy.optimize

G = 1534949.0910000005;
constant = 0.0244e7


class DeriveAGE:

    def __init__(self, interpolator, feh, mass, teff_value):
        """Set-up to use the given interpolator."""

        self.interp = interpolator
        self.feh = feh
        self.mass = mass
        self.teff_value = teff_value
        self.age_bound_check = False

    def teff_diff(self, age):

        """Return the effective temperature for the given stellar mass."""



        T = TeffK(self.interp('radius', self.mass, self.feh),
                  self.interp('lum', self.mass, self.feh)
                  )

        if age >T.min_age and age<T.max_age:

            return (T(age) - self.teff_value)

        else: return scipy.nan

    def possible_solution(self):

        """Finds 2 values of masses between which a solution is present"""

        age_array = scipy.linspace(0.01, 12, 100)
        teff_array_diff = scipy.empty(len(age_array))

        for a_index, a_value in enumerate(age_array):
            teff_array_diff[a_index] = self.teff_diff(a_value)
            print('for A = {}, dT = {}'.format(a_value,self.teff_diff(a_value)))

        #os.remove('age_test.txt')
        age_solutions = []
        for i in range(teff_array_diff.size - 1):
            x = teff_array_diff[i] * teff_array_diff[i + 1]
            with open('age_test.txt', 'a') as f:
                f.write(repr(i) + '\t' + repr(x) + '\n')
            if x < 0:
                age_solutions.append(age_array[i])
                age_solutions.append(age_array[i + 1])
                return age_solutions

        if i >=teff_array_diff.size - 2: self.age_bound_check = True

    def __call__(self):

        """Optimize the possible solution to calculate exact solution"""

        solution = 0
        age_solutions = self.possible_solution()
        print(age_solutions)

        if (self.age_bound_check) is False:
            solution = scipy.optimize.brentq(self.teff_diff, age_solutions[0],
                                             age_solutions[1])
            print("solution = ", solution)
            return solution
        else: return scipy.nan



if __name__=='__main__':


    system=sys.argv[1]



    with open('windemuth_stellar_raw.txt', 'r') as f:
        for lines in f:
            x=lines.split()
            if x[0]==system:
                feh=float(x[1])
                mass=float(x[7])
                cat_age=float(x[4])
                break

    with open('catalog_KIC.txt','r') as f:
        for lines in f:
            x=lines.split()
            if x[1]==system:
                teff=float(x[2])

    with open('age_test.txt','w') as f:
        f.write('i'+'\t'+'x'+'\n')

    serialized_dir = poet_path +  "/stellar_evolution_interpolators"
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    eccentricity_path=os.path.join(poet_path,'eccentricity_expansion_coef.txt').encode('ascii')

    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        eccentricity_path
    )


    age=DeriveAGE(interpolator,
                  feh,
                  mass,
                  teff)


    a=age()
    print('feh = ',feh)
    print('mass = ',mass)
    print('cat_age = ',cat_age)
    print('teff = ',teff)
    print('Derived Age = ',a)
