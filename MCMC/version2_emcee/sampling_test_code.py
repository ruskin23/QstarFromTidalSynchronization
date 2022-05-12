from os import system
from pickle import NONE
import numpy
import random
from numpy.core.fromnumeric import shape
import matplotlib.pyplot as plt
from numpy.core.function_base import linspace
import scipy
from scipy import special
import corner
import sys

import scipy.integrate as integrate

def get_cdf():

    M_array=priors[0]
    h_M=h[0]
    M_samples=all_chains[0]
    M_cdf=0
    for s in M_samples:
        alpha_M=(M_array-s)/h_M
        M_cdf=M_cdf+erf_fun(alpha_M)
    M_cdf=M_cdf/len(M_samples)

    Q_array=priors[1]
    h_Q=h[1]
    Q_samples=all_chains[1]
    Q_cdf=0
    for s in Q_samples:
        alpha_Q=(Q_array-s)/h_Q
        Q_cdf=M_cdf+erf_fun(alpha_Q)
    Q_cdf=Q_cdf/len(Q_samples)

    age_array=priors[3]
    h_age=h[3]
    age_samples=all_chains[3]
    age_cdf=0
    for s in age_samples:
        alpha_age=(M_array-s)/h_age
        age_cdf=age_cdf+erf_fun(alpha_age)
    age_cdf=age_cdf/len(age_samples)

    # z_pdf=get_pdf(z_array,z_samples,h_z)
    # plt.plot(z_array,z_pdf)
    # plt.show()
    # z_cdf=get_cdf(z_array,z_samples,h_z)

    # for h_z in [0.0002]:#,0.005,0.01,0.02,0.05,0.08,0.1]:
    #     z_array=priors[2]
    #     # h_M=h[0]
    #     z_samples=all_chains[2]
    #     z_pdf=get_pdf(z_array,z_samples,h_z)
    #     plt.plot(z_array,z_pdf)
    # plt.show()
   


def erf_fun(x):
    return 0.5*(1+special.erf(x/numpy.sqrt(2)))




def kernel_gauss(x,x_i,h):

    y=(x-x_i)/h
    return (1/(numpy.sqrt(2*numpy.pi)))*numpy.exp(-0.5*(y**2)) 

def kde(value,samples,bandwidth):

    n=len(samples)
    h=bandwidth

    f=0
    for i in range(n):
        f=f+kernel_gauss(value,samples[i],h)
    return f/(n*h)

def eccentricity_cdf(e_samples):

    # H=1.06*numpy.std(e_samples)*pow(len(e_samples),-0.2)
    H=0.01
    e_array=numpy.linspace(0,1,1000)

    e_cdf=0
    for s in e_samples:
        alpha_e=(e_array-s)/H
        e_cdf=e_cdf+erf_fun(alpha_e)
    e_cdf=e_cdf/len(e_samples)

    return e_cdf

    # plt.plot(e_array,e_cdf)
    # plt.show()

def eccentricity_pdf(e_samples):

        
    e_array=numpy.linspace(0,1,100)

    eccentricity_cdf(e_samples)

    print(numpy.std(e_samples))
    H=1.06*numpy.std(e_samples)*pow(len(e_samples),-0.2)
    print(H)
    
    F=kde(e_array,e_samples,H)
    F/=F.max()
    plt.plot(e_array,F,label=f'{H}')
    plt.legend()
    plt.show()

def eccentricity_sampling(rnd_num,e_samples):

    e_cdf=eccentricity_cdf(e_samples)
    e_array=numpy.linspace(0,1,10000)
    e_value=e_array[numpy.argmin(abs(e_cdf-rnd_num))]

    return e_value


if __name__=='__main__':

    system=sys.argv[1]

    A=numpy.load('/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/version2_emcee/catalog/samples/chains/'+system+'.npz')
    all_chains=numpy.transpose(A['thinned_chain'])

    esinw_samples=all_chains[8]
    ecosw_samples=all_chains[9]

    e_samples=[]
    for i in range(len(esinw_samples)):
        f=esinw_samples[i]**2 + ecosw_samples[i]**2
        f=numpy.sqrt(f)
        e_samples.append(f)
    e_samples=numpy.array(e_samples)
    

    sampling_parameters=['Msum','Q','z','age']

    priors=[numpy.linspace(0.2, 2.5, 5000),numpy.linspace(0.0085, 1, 5000),
    numpy.linspace(0.001, 0.06, 5000),10**(numpy.linspace(6, 10.1, 5000))/1e9]

    h=[0.01,0.001,0.0001,0.01]

    values=[]


    z_array=priors[2]
    h_z=h[2]
    z_samples=all_chains[2]
    z_cdf=0
    for s in z_samples:
        alpha_z=(z_array-s)/h_z
        z_cdf=z_cdf+erf_fun(alpha_z)
    z_cdf=z_cdf/len(z_samples)


    M_array=priors[0]
    h_M=h[0]
    M_samples=all_chains[0]

    Q_array=priors[1]
    h_Q=h[1]
    Q_samples=all_chains[1]


    age_array=priors[3]
    h_age=h[3]
    age_samples=all_chains[3]



    for i in range(5000):


        rdm_num=[random.uniform(0,1),random.uniform(0,1),random.uniform(0,1),random.uniform(0,1)]

        F_z=z_cdf
        z_value=z_array[numpy.argmin(abs(F_z-rdm_num[2]))]
        values.append(z_value)
        Z=numpy.exp(-0.5*(((z_value-z_samples)/h_z)**2))


        alpha_M=(M_array[None,:]-M_samples[:,None])/h_M

        F_M=numpy.dot(Z,erf_fun(alpha_M))
        F_M=F_M/F_M.max()
        M_value=M_array[numpy.argmin(abs(F_M-rdm_num[0]))]
        values.append(M_value)
        M=numpy.exp(-0.5*(((M_value-M_samples)/h_M)**2))

        alpha_Q=(Q_array[None,:]-Q_samples[:,None])/h_Q

        F_Q=numpy.dot(M*Z,erf_fun(alpha_Q))
        F_Q/=F_Q.max()
        Q_value=Q_array[numpy.argmin(abs(F_Q-rdm_num[1]))]
        values.append(Q_value)
        Q=numpy.exp(-0.5*(((Q_value-Q_samples)/h_Q)**2))

        alpha_age=(age_array[None,:]-age_samples[:,None])/h_age

        F_age=numpy.dot(Z*Q*M,erf_fun(alpha_age))
        F_age/=F_age.max()
        age_value=age_array[numpy.argmin(abs(F_age-rdm_num[3]))]
        values.append(age_value)

        e_rnd=random.uniform(0,1)

        values.append(eccentricity_sampling(e_rnd,e_samples))

    values=numpy.array(numpy.vstack(numpy.transpose(numpy.array(values))))
    d=values.reshape([len(values)//5,5])
    # numpy.savetxt('text.txt',d)
    figure=corner.corner(d)
    # plt.show()
    plt.savefig(f'corner_{system}.png')

