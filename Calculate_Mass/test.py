
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



class find_parameters():


    def teff_diff(self,mass,age):
        """Return the effective temperature for the given stellar mass."""

        T = TeffK(self.interp('radius', mass, self.feh),
                  self.interp('lum', mass, self.feh)
                  )

        try : return (T(age) - self.teff_value)
        except:
            #print('age error for mass: ', mass)
            return scipy.nan

    def find_range(self):

        T=scipy.empty([100,100])
        for i in range(100):
            for j in range(100):
                T[i,j]=self.teff_diff(self.mass_array[i],self.age_array[j])


        print(T[43])
        print(T[44])



    def __init__(self,
                 interpolator,
                 feh,
                 teff_value
                 ):

        self.feh=feh
        self.interp=interpolator
        self.teff_value=teff_value
        self.mass_array=numpy.linspace(0.4,1.2,100)
        self.age_array=numpy.linspace(1.1,5.2,100)




if __name__ == '__main__':

    serialized_dir = "/home/kpenev/projects/git/poet/stellar_evolution_interpolators"
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    teff=5580
    feh=-0.750

    T = find_parameters(interpolator,feh,teff)

    T.find_range()










