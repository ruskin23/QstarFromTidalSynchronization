#All Test Models between phi and logQ


import sys
sys.path.append('../../')
from utils import cummulative_distribution
from utils import Norm
from utils import multivariate_gaussian

import numpy



class TestModels:


    def prior_model(self,
                    parameter):

        p=[]

        if parameter=='phi':bounds=self.phi_bounds
        if parameter=='logQ':bounds=self.logQ_bounds
        for i in range(self.sample_space):
            p=numpy.append(m,random.uniform(bounds[0],bounds[1]))
        return cummulative_distribution(p)()


    """
    def gaussian_phi_model(self,
                           parameter,
                           phi_key):

        p=[]
        parameter_values=[]
        if parameter='phi':
            with open(self.sample_file,'r') as f:
                next(f)
                for lines in f:
                    m_i=float(x[3])
                    s_vector=[]
                    parameter_values=numpy.append(parameter_values,float(x[self.phi_keys.index(phi_key)]))
                    for key in self.phi_keys:
                        s_vector=numpy.append(s_vector,float(x[self.phi_keys.index(key)]))
                    p=numpy.append(m,m_i*multivariate_gaussian.multi*(s_vector,self.phi_mean_vector,s_vector,self.phi_mean_vector,self.phi_sigminv_matrix))
            cdf=cummulative_distribution(parameter_values,p)
            return cdf()
    """

    def gaussian_phi_model(self,
                           phi_key):


        phi_mean=self.phi_mean_vector[self.phi_keys.index(phi_key)]
        phi_sigma=self.phi_sigma_vector[self.phi_keys.index(phi_key)]

        phi_values=[]
        P=[]

        with open(self.samples_file,'r') as f:
            next(f)
            for lines in f:
                x=lines.split()
                phi_values=numpy.append(phi_values,float(x[self.phi_keys.index(phi_key)]))
                P=numpy.append(P,Norm.N(phi_value,loc=phi_mean,sigma=phi_sigma))

        cdf=cummulative_distribution(phi_values,P)
        return cdf()


    def gaussian_phi_phi_model(self,
                               phi,
                               phi_key1,
                               phi_key2):

        rho=self.rho[self.phi_keys.index(phi_key1),self.phi_keys.index(phi_key2)]
        phi1_mean=self.phi_mean_vector[self.phi_keys.index(phi_key1)]
        phi2_mean=self.phi_mean_vector[self.phi_keys.index(phi_key2)]
        phi1_sigma=self.phi_sigma_vector[self.phi_keys.index(phi_key1)]
        phi2_sigma=self.phi_sigma_vector[self.phi_keys.index(phi_key2)]


        phi_values=[]
        with open(self.samples_file,'r') as f:
            next(f)
            for lines in f:
                x=lines.split()
                phi_values=numpy.append(phi_values,float(x[self.phi_keys.index(phi)]))

        P=numpy.zeros(len(phi_values))
        with open(self.samples_file,'r') as f:
            next(f)
            for lines in f:
                x=lines.split()
                phi2=float(x[self.phi_keys.index(phi_key2)])
                P=P+Norm.N(phi_values,loc=phi1_mean,sigma=phi1_sigma)*Norm.N(phi2,loc=phi2_mean,sigma=phi2_sigma)*numpy.exp(-(phi_values-phi1_mean)*(phi2-phi2_mean)/(rho*phi1_sigma*phi2_sigma))

        cdf=cummulative_distribution(phi_values,P)
        return cdf()



    def integrand(self,
                  phi,
                  logQ,
                  logQ_mean,
                  logQ_sigma,
                  phi_mean):

        return Norm.N(logQ,loc=logQ_mean,scale=logQ_sigma)*numpy.exp(-(logQ-logQ_mean)*(phi-phi_mean)/(self.rho*logQ_mean*phi_mean))

    def Porbablity(self,
                   phi,
                   logQ_mean,
                   logQ_sigma,
                   phi_mean,
                   phi_sigma,
                   bounds):

        I=lambda logQ: self.integrand(phi,logQ,logQ_mean,logQ_sigma,phi_mean)
        Integral = integrate.quad(I,bounds[0],bounds[1])[0]

        return Norm.N(phi,loc=phi_mean,scale=phi_sigma)


    def correlated_phi_theta(self,
                             parameter,
                             phi_key):


        phi_mean=self.phi_mean_vector[self.phi_keys.index(phi_key)]
        phi_sigma=self.phi_sigma_vector[self.phi_keys.index(phi_key)]

        phi_logQ_sigma=self.rho[len(self.phi_keys),self.phi_keys.index(phi_key)]*phi_sigma*self.logQ_mean

        bounds=self.logQ_bounds

        if parameter=='logQ':
            logQ_samples=numpy.linspace(bound[0],bound[1],1000)
            P=numpy.zeros(len((logQ_samples)))
            with open(self.sample_file,'r') as f:
                next(f)
                for lines in f:
                    x=lines.split()
                    phi_value=float(x[self.phi_keys.index(phi_key)])
                    P=P+(Norm.N(phi_value,loc=phi_mean,sigma=phi_sigma)*Norm.N(logQ_samples,loc=self.logQ_mean,sigma=self.logQ_sigma)*numpy.exp(-0.5*(phi_value-phi_mean)*(logQ_samples-self.logQ_mean)/phi_logQ_sigma))
            cdf=cummulative_distribution(logQ_samples,P)()
            return cdf

        if parameter=='phi':

            P=[]
            phi_values=[]
            with open(self.sample_file,'r') as f:
                next(f)
                for lines in f:
                    x=lines.split()
                    phi_value=float(x[self.phi_keys.index(phi_key)])
                    phi_values=numpy.append(phi_values,phi_value)
                    P=numpy.append(P,self.Probability(phi_value,self.logQ_mean,self.logQ_sigma,phi_mean,phi_sigma,bounds))

            cdf=cummulative_distribution(phi_values,P)()

            return cdf



    def __init__(self,
                 samples_file,
                 logQ_mean,
                 logQ_sigma,
                 logQ_bounds,
                 rho,
                 phi_keys,
                 phi_mean_vector,
                 phi_sigma_vector,
                 phi_bounds=None):

        self.samples_file=samples_file
        self.logQ_mean=logQ_mean
        self.logQ_sigma=logQ_sigma
        self.logQ_bounds=logQ_bounds

        self.phi_keys=phi_keys,
        self.phi_mean_vector=phi_mean_vector
        self.phi_sigma_vector=self.phi_sigma_vector

        d=numpy.diag_indices(3)
        self.phi_sigmainv_matrix=numpy.zeros(len(phi_keys),len(phi_keys))
        self.phi_sigmainv_matrix[d]=1/(phi_mean_vector**2)

        self.rho=rho
