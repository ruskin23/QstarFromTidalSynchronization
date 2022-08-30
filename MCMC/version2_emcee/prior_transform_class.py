import numpy
from utils import erf_fun
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


    # def paramter_evaluate(self,uniform_values):

    #     Z_cdf=0.0
    #     for s in self.Z_samples:
    #         z=(self.Z_prior-s)/self.h_Z
    #         Z_cdf+=erf_fun(z)
    #     Z_cdf/=len(self.Z_samples)
    #     Z_cdf_interp = interp1d(Z_cdf,self.Z_prior)

    #     if uniform_values[0]<min(Z_cdf):Z_value=min(self.Z_samples)
    #     elif uniform_values[0]>max(Z_cdf):Z_value=max(self.Z_samples)
    #     else:Z_value=Z_cdf_interp(uniform_values[0])
    #     Z_p=numpy.exp(-0.5*(((Z_value-self.Z_samples)/self.h_Z)**2))

    #     q=(self.Q_prior[None,:]-self.Q_samples[:,None])/self.h_Q
    #     erfq=erf_fun(q)
    #     F_Q=numpy.dot(Z_p,erfq)
    #     F_Q/=F_Q.max()
    #     Q_value=self.Q_prior[numpy.argmin(abs(F_Q-uniform_values[1]))]
    #     Q_p=numpy.exp(-0.5*(((Q_value-self.Q_samples)/self.h_Q)**2))

    #     m=(self.M_prior[None,:]-self.M_samples[:,None])/self.h_M
    #     erfm=erf_fun(m)
    #     F_M=numpy.dot(Z_p*Q_p,erfm)
    #     F_M/=F_M.max()
    #     M_value=self.M_prior[numpy.argmin(abs(F_M-uniform_values[2]))]
    #     M_p=numpy.exp(-0.5*(((M_value-self.M_samples)/self.h_M)**2))

    #     t=(self.t_prior[None,:]-self.t_samples[:,None])/self.h_t
    #     erft=erf_fun(t)
    #     F_t=numpy.dot(Z_p*Q_p*M_p,erft)
    #     F_t/=F_t.max()
    #     t_value=self.t_prior[numpy.argmin(abs(F_t-uniform_values[3]))]
    #     t_p=numpy.exp(-0.5*(((t_value-self.t_samples)/self.h_t)**2))

    #     e=(self.e_prior[None,:]-self.e_samples[:,None])/self.h_e
    #     erfe=erf_fun(e)
    #     F_e=numpy.dot(Z_p*Q_p*M_p*t_p,erfe)
    #     F_e/=F_e.max()
    #     e_value=self.e_prior[numpy.argmin(abs(F_e-uniform_values[4]))]

    #     return [Z_value,M_value,Q_value,t_value,e_value] 
