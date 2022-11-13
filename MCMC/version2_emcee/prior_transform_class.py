import numpy
import logging
import pickle

from utils import DiscreteSampling

from pathlib import Path
from directories import directories

home_dir=str(Path.home())
path=directories(home_dir)


import logging

_logger=logging.getLogger(__name__)

class prior_transform:

    def __init__(self,system_num):

        self.system_num=system_num

        system_chains=numpy.load(path.current_directory+'/catalog/samples/chains/'+system_num+'.npz')
        all_chains=numpy.transpose(system_chains['thinned_chain'])

        M_samples=all_chains[0]
        Q_samples=all_chains[1]
        Z_samples=all_chains[2]
        t_samples=all_chains[3]

        esinw_samples=all_chains[8]
        ecosw_samples=all_chains[9]
        
        e_samples=numpy.sqrt(ecosw_samples**2 + esinw_samples**2)

        with open(path.current_directory+'/bandwidth.pickle','rb') as f:
            bw=pickle.load(f)
        bw=bw[system_num]['ISJ']
        _logger.info('Bandwidths: {!r}'.format(bw))
        with open(path.current_directory+'/e_bandwidth.pickle','rb') as f:
            bw_e=pickle.load(f)
        _logger.info('Eccentricity Bandwidth: {!r}'.format(bw_e[system_num]))

        h_M=bw[0]
        h_Q=bw[1]
        h_Z=bw[2]
        h_t=bw[3]
        h_e=bw_e[system_num]

        self._quantities=['Z','M','Q','t','e']

        self._samples={'M':M_samples,
                       'Q':Q_samples,
                       'Z':Z_samples,
                       't':t_samples,
                       'e':e_samples}

        self._bandwidths={'M':h_M,
                          'Q':h_Q,
                          'Z':h_Z,
                          't':h_t,
                          'e':h_e}


    def get_samples(self,unit_cube):

        weights=numpy.ones(len(self._samples[self._quantities[0]]))
        unit_cube_iter=iter(unit_cube)

        sampled_value=[]

        for Q in self._quantities:

            DS=DiscreteSampling(self._samples[Q],self._bandwidths[Q])

            DS.set_weights(weights)

            sampled=DS.ppf(next(unit_cube_iter))

            weights*=DS.get_norm(sampled)

            sampled_value.append(sampled)
        
        return sampled_value



