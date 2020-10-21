
import sys

from pathlib import Path
home_dir=str(Path.home())

if home_dir=='/home/rxp163130':
    poet_path=home_dir+'/poet'
if home_dir=='/home/ruskin':
    poet_path=home_dir+'/projects/poet'
if home_dir=='/home1/06850/rpatel23':
    work_dir='/work/06850/rpatel23/stampede2'
    poet_path=work_dir+'/poet'

sys.path.append(poet_path+'/PythonPackage')
sys.path.append(poet_path+'/scripts')

import sys

from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from orbital_evolution.binary import Binary
from orbital_evolution.transformations import phase_lag
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet
from stellar_evolution.manager import StellarEvolutionManager
from stellar_evolution.derived_stellar_quantities import\
    TeffK,\
    LogGCGS,\
    RhoCGS
