import numpy
from utils import erf_fun
import logging

from scipy.interpolate import interp1d

from pathlib import Path
from directories import directories

home_dir=str(Path.home())
path=directories(home_dir)

class prior_transform:

    def __init__(self,system_num,sampling_util_dict):

        self.system_num=system_num

        for item,value in sampling_util_dict.items():
            setattr(self,item,value)

        # system_chains=numpy.load(path.current_directory+'/catalog/samples/chains/'+system_num+'.npz')
        # all_chains=numpy.transpose(system_chains['thinned_chain'])

        # self.M_samples=all_chains[0]
        # self.Q_samples=all_chains[1]
        # self.Z_samples=all_chains[2]
        # self.t_samples=all_chains[3]

        # ecosw_samples=all_chains[8]
        # esinw_samples=all_chains[9]
        # self.e_samples=numpy.sqrt(ecosw_samples**2 + esinw_samples**2)

        # self.M_prior=numpy.linspace(min(self.M_samples), max(self.M_samples), 10000)
        # self.Q_prior=numpy.linspace(min(self.Q_samples), max(self.Q_samples), 10000)
        # self.Z_prior=numpy.linspace(min(self.Z_samples), max(self.Z_samples),10000)
        # self.t_prior=numpy.linspace(min(self.t_samples), max(self.t_samples), 10000)
        # self.e_prior=numpy.linspace(min(self.e_samples),max(self.e_samples),10000)

        # self.bandwidths=bandwidths

        # self.h_M=bw_tuple[0]
        # self.h_Q=bw_tuple[1]
        # self.h_Z=bw_tuple[2]
        # self.h_t=bw_tuple[3]
        # self.h_e=bw_tuple[4]


    def paramter_evaluate(self,uniform_values):

        Z_cdf=0.0
        for s in self.Z_samples:
            z=(self.Z_prior-s)/self.h_Z
            Z_cdf+=erf_fun(z)
        Z_cdf/=len(self.Z_samples)
        Z_cdf_interp = interp1d(Z_cdf,self.Z_prior)

        if uniform_values[0]<min(Z_cdf):Z_value=min(self.Z_samples)
        elif uniform_values[0]>max(Z_cdf):Z_value=max(self.Z_samples)
        else:Z_value=Z_cdf_interp(uniform_values[0])
        Z_p=numpy.exp(-0.5*(((Z_value-self.Z_samples)/self.h_Z)**2))

        # q=(self.Q_prior[None,:]-self.Q_samples[:,None])/self.h_Q
        # erfq=erf_fun(q)
        F_Q=numpy.dot(Z_p,self.erfq)
        F_Q/=F_Q.max()
        Q_value=self.Q_prior[numpy.argmin(abs(F_Q-uniform_values[1]))]
        Q_p=numpy.exp(-0.5*(((Q_value-self.Q_samples)/self.h_Q)**2))

        # m=(self.M_prior[None,:]-self.M_samples[:,None])/self.h_M
        # erfm=erf_fun(m)
        F_M=numpy.dot(Z_p*Q_p,self.erfm)
        F_M/=F_M.max()
        M_value=self.M_prior[numpy.argmin(abs(F_M-uniform_values[2]))]
        M_p=numpy.exp(-0.5*(((M_value-self.M_samples)/self.h_M)**2))

        # t=(self.t_prior[None,:]-self.t_samples[:,None])/self.h_t
        # erft=erf_fun(t)
        F_t=numpy.dot(Z_p*Q_p*M_p,self.erft)
        F_t/=F_t.max()
        t_value=self.t_prior[numpy.argmin(abs(F_t-uniform_values[3]))]
        t_p=numpy.exp(-0.5*(((t_value-self.t_samples)/self.h_t)**2))

        # e=(self.e_prior[None,:]-self.e_samples[:,None])/h_e
        # erfe=erf_fun(e)
        F_e=numpy.dot(Z_p*Q_p*M_p*t_p,self.erfe)
        F_e/=F_e.max()
        e_value=self.e_prior[numpy.argmin(abs(F_e-uniform_values[4]))]

        return [Z_value,M_value,Q_value,t_value,e_value] 
