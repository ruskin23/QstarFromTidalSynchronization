#All Test Models between phi and logQ

import sys
sys.path.append('../../')
from utils import cummulative_distribution
from utils import Norm
from utils import multivariate_gaussian

import numpy
from scipy import integrate
from scipy import special
import random

class TestModels:


    def prior_model(self,
                    phi_key,
                    parameter):

        p=[]

        if parameter=='phi':bounds=numpy.array([self.sampling_parameters[phi_key]['min'],self.sampling_parameters[phi_key]['max']])
        if parameter=='logQ':bounds=self.logQ_bounds
        for i in range(1000):
            p=numpy.append(p,random.uniform(bounds[0],bounds[1]))
        return cummulative_distribution(p)()


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
                m_i=float(x[3])
                phi_value=float(x[self.phi_keys.index(phi_key)])
                phi_values=numpy.append(phi_values,phi_value)
                P=numpy.append(P,m_i*Norm().N(phi_value,loc=phi_mean,sigma=phi_sigma))

        print(numpy.mean(phi_values))
        cdf=cummulative_distribution(phi_values,P)
        return cdf()


    def gaussian_phi_theta_model(self,
                                 parameter,
                                 phi_key):


        phi_mean=self.phi_mean_vector[self.phi_keys.index(phi_key)]
        phi_sigma=self.phi_sigma_vector[self.phi_keys.index(phi_key)]

        print(phi_mean)
        print(phi_sigma)

        P=[]

        if parameter in self.phi_keys:return self.gaussian_phi_model(phi_key)

        if parameter=='logQ':
            print('MEAN  = ',self.logQ_mean)
            logQ_values=numpy.linspace(self.logQ_bounds[0],self.logQ_bounds[1],1000)
            P=Norm().N(logQ_values,loc=self.logQ_mean,sigma=self.logQ_sigma)

            cdf=cummulative_distribution(logQ_values,P)
            return cdf()


    def Probability(self,
                   vN,
                   x_mean,
                   x_sigma,
                   vN_mean,
                   vN_sigma,
                   bounds):


        I=lambda x:Norm().N(x,loc=x_mean,sigma=x_sigma)*numpy.exp(-(x-x_mean)*(vN-vN_mean)/(self.rho*x_sigma*vN_sigma))
        Integral = integrate.quad(I,bounds[0],bounds[1])[0]
        return Norm().N(vN,loc=vN_mean,sigma=vN_sigma)*Integral


    def correlated_phi_theta(self,
                             parameter,
                             phi_key):

        phi_mean=self.phi_mean_vector[self.phi_keys.index(phi_key)]
        phi_sigma=self.phi_sigma_vector[self.phi_keys.index(phi_key)]
        phi_logQ_sigma=self.rho*phi_sigma*self.logQ_sigma

        bounds=self.logQ_bounds

        if parameter=='logQ':
            logQ_samples=numpy.linspace(bounds[0],bounds[1],1000)
            P=numpy.zeros(len((logQ_samples)))
            with open(self.samples_file,'r') as f:
                next(f)
                for lines in f:
                    x=lines.split()
                    phi_value=float(x[self.phi_keys.index(phi_key)])
                    arg1=(phi_value-phi_sigma)/phi_sigma
                    arg2=(logQ_samples-self.logQ_mean)/self.logQ_sigma
                    arg3=2*self.rho*arg1*arg2
                    P=P+numpy.exp(-0.5*(1/(1-self.rho**2))*(arg1**2 + arg2**2 - arg3))
                    #P=P+(Norm().N(phi_value,loc=phi_mean,sigma=phi_sigma)*Norm().N(logQ_samples,loc=self.logQ_mean,sigma=self.logQ_sigma)*numpy.exp(-(phi_value-phi_mean)*(logQ_samples-self.logQ_mean)/phi_logQ_sigma))
            cdf=cummulative_distribution(logQ_samples,P)()
            return cdf

        if parameter in self.phi_keys:

            P=[]
            phi_values=[]
            N=0
            with open(self.samples_file,'r') as f:
                next(f)
                for lines in f:
                    x=lines.split()
                    m_i=float(x[3])
                    phi_value=float(x[self.phi_keys.index(phi_key)])
                    phi_values=numpy.append(phi_values,phi_value)
                    b=(phi_value-phi_mean)/phi_sigma
                    x=(self.logQ_bounds-self.logQ_mean)/self.logQ_sigma
                    y=x-self.rho*b
                    z=y/(1-self.rho**2)
                    p=m_i*numpy.exp(-0.5*b*b)*(special.erf(z[1])-special.erf(z[0]))
                    #p=m_i*self.Probability(phi_value,self.logQ_mean,self.logQ_sigma,phi_mean,phi_sigma,bounds)
                    P=numpy.append(P,p)
                    N=N+p
            P=P/N
            cdf=cummulative_distribution(phi_values,P)()

            return cdf



    def gaussian_theta_theta_model(self,
                                   theta_key,
                                   parameter):


        if parameter=='logQ':

            logQ_samples=numpy.linspace(self.logQ_bounds[0],self.logQ_bounds[1],1000)
            p=numpy.zeros(1000)
            for i,q in enumerate(logQ_samples):
                p[i]=Norm().N(q,loc=self.logQ_mean,sigma=self.logQ_sigma)
            cdf=cummulative_distribution(logQ_samples,p)()
            return cdf

        else:
            theta_mean=self.sampling_parameters[theta_key]['value']
            theta_sigma=self.model_width[theta_key]
            if self.sampling_parameters[theta_key]['dist']=='Uniform':
                bounds=[self.sampling_parameters[theta_key]['min'],self.sampling_parameters[theta_key]['max']]
            else:
                #bounds=[theta_mean-5*theta_sigma,theta_mean+5*theta_sigma]
                bounds =[0,1]
            theta_samples=numpy.linspace(bounds[0],bounds[1],1000)
            p=numpy.zeros(1000)
            for i,q in enumerate(theta_samples):
                if q>0:
                    p[i]=Norm().N(q,loc=theta_mean,sigma=theta_sigma)*Norm().N(q,loc=theta_mean,sigma=self.sampling_parameters[theta_key]['sigma'])
                else:p[i]=0
            cdf=cummulative_distribution(theta_samples,p)()
            return cdf


    def correlated_theta_theta(self,
                               theta_key,
                               parameter):

        theta_mean=self.sampling_parameters[theta_key]['value']
        theta_sigma=self.model_width[theta_key]

        if parameter=='logQ':
            logQ_samples=numpy.linspace(self.logQ_bounds[0],self.logQ_bounds[1],1000)
            p=numpy.zeros(1000)
            if self.sampling_parameters[theta_key]['dist']=='Uniform':
                bounds=[self.sampling_parameters[theta_key]['min'],self.sampling_parameters[theta_key]['max']]
                for i,q in enumerate(logQ_samples):
                    p[i]=self.Probability(q,theta_mean,theta_sigma,self.logQ_mean,self.logQ_sigma,bounds)
                cdf=cummulative_distribution(logQ_samples,p)()
                return cdf
            else:
                for i,q in enumerate(logQ_samples):
                    b=(q-self.logQ_mean)/self.logQ_sigma
                    theta_lim=numpy.array([0,1])
                    x=(theta_lim-theta_mean)/theta_sigma
                    y=x-self.rho*b
                    z=y/(1-self.rho**2)
                    p[i]=numpy.exp(-b**2/2)*(special.erf(z[1])-special.erf(z[0]))
                    #arg=(q-self.logQ_mean)/self.logQ_sigma
                    #l=(arg/self.rho) - (theta_mean/theta_sigma)
                    #h=(arg/self.rho) + ((1-theta_mean)/theta_sigma)
                    #Integral=numpy.sqrt(2*numpy.pi)*(special.erf(h)-special.erf(l))
                    #p[i]=numpy.exp(-0.5*arg*arg*(1-(1/self.rho)**2))*theta_sigma*Integral
                    #print('q={},q_0={},sigma={},qrg={},l={},h={},i={}'.format(q,self.logQ_mean,self.logQ_sigma,arg,l,h,Integral))
                    #I=lambda x:numpy.exp(-0.5*x*x)
                    #Integral = integrate.quad(I,l,h)[0]
                    #p[i]=self.Probability(q,theta_mean,theta_sigma,self.logQ_mean,self.logQ_sigma,[0,0.5])
                    #p[i]=Norm().N(q,loc=self.logQ_mean,sigma=self.logQ_sigma)
                cdf=cummulative_distribution(logQ_samples,p)()
                print('Q_mean = ',self.logQ_mean)
                return cdf


        else:

            if self.sampling_parameters[theta_key]['dist']=='Uniform':
                bounds=[self.sampling_parameters[theta_key]['min'],self.sampling_parameters[theta_key]['max']]
            else:
                bounds =[-10*theta_mean,10*theta_mean]
            theta_samples=numpy.linspace(bounds[0],bounds[1],1000)
            p=numpy.zeros(1000)
            for i,q in enumerate(theta_samples):
                if q>0:
                    b=(q-theta_mean)/theta_sigma
                    theta_lim=numpy.array([5,12])
                    x=(theta_lim-self.logQ_mean)/self.logQ_sigma
                    y=x-self.rho*b
                    z=y/(1-self.rho**2)
                    p[i]=numpy.exp(-b**2/2)*(special.erf(z[1])-special.erf(z[0]))*Norm().N(q,loc=theta_mean,sigma=self.sampling_parameters[theta_key]['sigma'])
                    #p[i]=self.Probability(q,self.logQ_mean,self.logQ_sigma,theta_mean,theta_sigma,self.logQ_bounds)*Norm().N(q,loc=theta_mean,sigma=self.sampling_parameters[theta_key]['sigma'])
                else:p[i]=0
            cdf=cummulative_distribution(theta_samples,p)()
            return cdf





    def __init__(self,
                 samples_file,
                 sampling_parameters,
                 model_width,
                 rho):



        self.samples_file=samples_file
        self.sampling_parameters=sampling_parameters
        self.model_width=model_width
        self.rho=rho

        self.logQ_mean=self.sampling_parameters['logQ']['value']
        self.logQ_sigma=self.model_width['logQ']
        self.logQ_bounds=numpy.array([self.sampling_parameters['logQ']['min'],self.sampling_parameters['logQ']['max']])


        self.phi_keys=[]
        self.phi_mean_vector=[]
        self.phi_sigma_vector=[]
        for key,value in self.sampling_parameters.items():
            if value['dist']=='Samples':
                self.phi_keys.append(key)
                self.phi_mean_vector=numpy.append(self.phi_mean_vector,value['value'])
                self.phi_sigma_vector=numpy.append(self.phi_sigma_vector,self.model_width[key])


