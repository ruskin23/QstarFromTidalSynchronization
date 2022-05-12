import sys
import os
import random
import numpy
sys.path.append('/home/ruskin/projects/poet/PythonPackage')
sys.path.append('/home/ruskin/projects/poet/scripts')
sys.path.append('/home/ruskin/projects/QstarFromTidalSynchronization/binary_star_evolution/updated_evolution_code')

from stellar_evolution.library_interface import library

from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as \
    orbital_evolution_library
from orbital_evolution.transformations import phase_lag
from stellar_evolution.derived_stellar_quantities import\
    TeffK,\
    LogGCGS,\
    RhoCGS


from create_objects import BinaryObjects


if __name__=='__main__':


    serialized_dir =  "/home/ruskin/projects/poet/stellar_evolution_interpolators"
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    eccentricity_path=os.path.join('/home/ruskin/projects/poet','eccentricity_expansion_coef.txt').encode('ascii')

    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        eccentricity_path
    )

    with open('catalog/filtering/normandy_poet_catalog.txt','r') as f:
        next(f)
        for i,lines in enumerate(f):
            x=lines.split()
            KIC=x[1]
            Porb=float(x[2])
            e=float(x[3])
            z=float(x[4])
            age=float(x[5])
            m1=float(x[6])
            m2=float(x[7])
            feh=library.feh_from_z(z)

            T=10**(age-9)

            # print(f'{KIC}\t{m1}\t{m2}\t{T}\t{feh}')

            try:
                primary_radius=interpolator('radius',m1,feh)
                primary_lum=interpolator('lum',m1,feh)
                G=LogGCGS(m1,primary_radius)
                primary_logg=G(T)

                secondary_radius=interpolator('radius',m2,feh)
                secondary_lum=interpolator('lum',m2,feh)
                G=LogGCGS(m1,secondary_radius)
                secondary_logg=G(T)

                print(f'{KIC}\t{primary_logg}\t{secondary_logg}')
            except: pass







