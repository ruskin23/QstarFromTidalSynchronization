import numpy
import sys

from pathlib import Path
from directories import directories

home_dir=str(Path.home())
path=directories(home_dir)

sys.path.append(path.git_dir)

from utils import DiscreteSampling, _get_kernel_bandwidth

class discrete_stats:

    def __init__(self,samples):

        self._flat_posteriors = {}
        for name, values in samples:
            self._flat_posteriors[name] = {'samples': values,
                                            'bandwidth': _get_kernel_bandwidth(values)}
            assert(len(values) == self._number_samples)
            self._number_samples = len(values)

        self.weights = numpy.ones(self._number_samples )        

    def _assign_pdf(self, rv_vector):

        rv_vector_iter = iter(rv_vector)

        for params, dist in self._flat_posteriors.items():

            distrbution = DiscreteSampling(dist['samples'], dist['bandwidth'])
            distrbution.set_weights(self.weights)
            dist['pdf'] = self.weights*distrbution._pdf(next(rv_vector_iter))
            self.weights = distrbution.get_norm(next(rv_vector))

        







