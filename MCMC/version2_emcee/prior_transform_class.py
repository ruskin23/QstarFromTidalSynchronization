import numpy
from utils import erf_fun
import random
import corner
import matplotlib.pyplot as plt
import logging

_logger = logging.getLogger(__name__)

class prior_transform:

    def __init__(self,system_num):

        self.system_num=system_num

        system_chains=numpy.load('/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/version2_emcee/catalog/samples/chains/'+system_num+'.npz')
        all_chains=numpy.transpose(system_chains['thinned_chain'])

        self.M_samples=all_chains[0]
        self.Q_samples=all_chains[1]
        self.z_samples=all_chains[2]
        self.t_samples=all_chains[3]

        ecosw_samples=all_chains[8]
        esinw_samples=all_chains[9]
        self.e_samples=numpy.sqrt(ecosw_samples**2 + esinw_samples**2)


        self.M_prior=numpy.linspace(0.2, 2.5, 1000)
        self.Q_prior=numpy.linspace(0.0085, 1, 1000)
        self.z_prior=numpy.linspace(0.001, 0.06, 1000)
        self.t_prior=10**(numpy.linspace(6, 10.1, 1000))/1e9
        self.e_prior=numpy.linspace(0,1,1000)

        self.bandwidth=[0.01,0.001,0.0001,0.001]

        self.h_M=self.bandwidth[0]
        self.h_Q=self.bandwidth[1]
        self.h_z=self.bandwidth[2]
        self.h_t=self.bandwidth[3]
        self.h_e=0.01

    def paramter_evaluate(self,uniform_values):

        M_cdf=0.0
        for s in self. M_samples:
            m=(self.M_prior-s)/self.h_M
            M_cdf+=erf_fun(m)
        M_cdf/=len(self.M_samples)
        M_value=self.M_prior[numpy.argmin(abs(M_cdf-uniform_values[0]))]
        M_p=numpy.exp(-0.5*(((M_value-self.M_samples)/self.h_M)**2))

        q=(self.Q_prior[None,:]-self.Q_samples[:,None])/self.h_Q
        F_Q=numpy.dot(M_p,erf_fun(q))
        F_Q/=F_Q.max()
        Q_value=self.Q_prior[numpy.argmin(abs(F_Q-uniform_values[1]))]
        Q_p=numpy.exp(-0.5*(((Q_value-self.Q_samples)/self.h_Q)**2))

        z=(self.z_prior[None,:]-self.z_samples[:,None])/self.h_z
        F_z=numpy.dot(Q_p*M_p,erf_fun(z))
        F_z/=F_z.max()
        z_value=self.z_prior[numpy.argmin(abs(F_z-uniform_values[2]))]
        z_p=numpy.exp(-0.5*(((z_value-self.z_samples)/self.h_z)**2))

        t=(self.t_prior[None,:]-self.t_samples[:,None])/self.h_t
        F_t=numpy.dot(z_p*M_p*Q_p,erf_fun(t))
        F_t/=F_t.max()
        t_value=self.t_prior[numpy.argmin(abs(F_t-uniform_values[3]))]

        e_cdf=0
        for s in self.e_samples:
            alpha_e=(self.e_prior-s)/self.h_e
            e_cdf=e_cdf+erf_fun(alpha_e)
        e_cdf=e_cdf/len(self.e_samples)

        e_value=self.e_prior[numpy.argmin(abs(e_cdf-uniform_values[4]))]

        return [M_value,Q_value,z_value,t_value,e_value]


def corner_plot(system_number):

    prior=prior_transform(system_number)

    values=[]
    for i in range(1000):
        rnd_nums=[random.uniform(0,1),random.uniform(0,1),random.uniform(0,1),random.uniform(0,1),random.uniform(0,1)]

        values.append(prior.paramter_evaluate(rnd_nums))

    flat_values = [item for sublist in values for item in sublist]

    d=numpy.array(numpy.vstack(numpy.transpose(numpy.array(flat_values))))
    data_corner=d.reshape([len(d)//5,5])
    figure=corner.corner(data_corner)
    # plt.show()
    plt.savefig('test_corner.png')




if __name__ == '__main__':
    system_number='10198109'

    # corner_plot(system_number)

    rnd_nums=[random.uniform(0,1),random.uniform(0,1),random.uniform(0,1),random.uniform(0,1),random.uniform(0,1)]


    sampler=prior_transform(system_number)

    print(sampler.paramter_evaluate(rnd_nums))
        







