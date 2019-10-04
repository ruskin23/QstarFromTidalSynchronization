import numpy

import sys
import os
import os.path

from pathlib import Path
home_dir=str(Path.home())

if home_dir=='/home/rxp163130':
    poet_path=home_dir+'/poet/'
if home_dir=='/home/ruskin':
    poet_path=home_dir+'/projects/poet/'

sys.path.append(poet_path+'PythonPackage')


from stellar_evolution.manager import StellarEvolutionManager
from stellar_evolution.derived_stellar_quantities import\
    TeffK,\
    LogGCGS,\
    RhoCGS


import matplotlib.pyplot as plt

class Quantity_Evaluator:

    def __init__(self,
                 interpolator,
                 mass,
                 feh):

        self.interpolator=interpolator
        self.mass=mass
        self.feh=feh

        self.luminosity=self.interpolator('lum',self.mass,self.feh)
        self.radius=self.interpolator('radius',self.mass,self.feh)


    def calculate_quantity(self,
                           name,
                           age):



        if name=='teff':
            T=TeffK(self.radius,self.luminosity)
            try:
                return T(age)
            except:
                return numpy.nan

        elif name=='logg':
            G=LogGCGS(self.mass,self.radius)
            try:
                return G(age)
            except:
                return numpy.nan

        elif name=='rho':
            R=RhoCGS(self.mass,self.radius)
            return R(age)

        else:
            raise ValueError('quantity not found')



if __name__ == '__main__':

    serialized_dir = poet_path +  "stellar_evolution_interpolators"
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    N=100
    mass_array=numpy.linspace(0.5,1.2,N)
    age_array=numpy.linspace(00.01,12,N)

    M=[0.6,0.7,1.0]

    for m in M:
        x=Quantity_Evaluator(interpolator,m,0.0)
        T=x.calculate_quantity('teff',age_array)
        print('MASS = ',m)
        print(T)



