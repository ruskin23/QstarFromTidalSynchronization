
import sys
sys.path.append('/Users/ruskinpatel/Desktop/Research/poet/PythonPackage')
sys.path.append('/Users/ruskinpatel/Desktop/Research/poet/scripts')


from stellar_evolution.manager import StellarEvolutionManager


from mass_calculations import DerivePrimnaryMass
from mass_calculations import DeriveSecondaryMass

serialized_dir = "/Users/ruskinpatel/Desktop/Research/poet/stellar_evolution_interpolators"
manager = StellarEvolutionManager(serialized_dir)
interpolator = manager.get_interpolator_by_name('default')


Pmass = DerivePrimnaryMass(interpolator,0.1,5.0,5886)
x = Pmass()




