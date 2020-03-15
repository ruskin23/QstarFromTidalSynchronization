import numpy
import scipy
import sys
sys.path.append('../../')
from utils import multivariate_gaussian
import random

class UncorrelatedSampling:


    def stepping_function_normalization(self,
                                        phi_vector):

        P=0

        with open(self.samples_file,'r') as f:
            next(f)
            for lines in f:
                x=lines.split()
                m_i=float(x[3])
                s_vector=[]
                for key in self.phi_keys:
                    s_vector=numpy.append(s_vector,float(x[self.phi_keys.index(key)]))
                P=P+m_i*multivariate_gaussian().multi(s_vector,phi_vector,s_vector,phi_vector,self.phi_sigma_matrix)

        return P

    def get_phi_samples(self):

        N=self.stepping_function_normalization(self.phi_mean)
        mean_vector=self.phi_mean
        sigma_vector=self.phi_sigma

        samples=[]
        modified_multiplicity=[]

        P=0
        with open(self.samples_file,'r') as f:
            next(f)
            for lines in f:
                x=lines.split()
                m_i=float(x[3])
                s_vector=[]
                for key in self.phi_keys:
                    s_vector=numpy.append(s_vector,float(x[self.phi_keys.index(key)]))
                samples.append(s_vector)
                modified_multiplicity.append(m_i*multivariate_gaussian().multi(s_vector,mean_vector,s_vector,mean_vector,self.phi_sigma_matrix))
                #arg_vector=((s_vector-mean_vector)/sigma_vector)**2
                #arg=0
                #for a in arg_vector:
                #    arg=arg+a
                #samples.append(s_vector)
                #modified_multiplicity.append(m_i*numpy.exp(0.5*arg))

        modified_multiplicity=modified_multiplicity/N
        U=random.uniform(0,1)
        for mm,s in zip(modified_multiplicity,samples):
            if mm>U:
                for key in self.phi_keys:
                    self.proposed[key]=s[self.phi_keys.index(key)]
                print('Proposed im sampling class = ',self.proposed)
                break
            else:U=U-mm

    def get_theta_samples(self):

        for key,value in self.sampling_parameters.items():
            if value['dist']=='Normal' or value['dist']=='Uniform':
                self.proposed[key]=scipy.stats.norm.rvs(loc=self.current_parameters[key],scale=value['step'])

    def __call__(self):

        self.get_theta_samples()
        self.get_phi_samples()
        return self.proposed


    def __init__(self,
                 system,
                 sampling_parameters,
                 parameters):

        self.system=system
        self.sampling_parameters=sampling_parameters
        self.current_parameters=parameters

        self.samples_file='../../../mcmc_mass_age/samples/updated_samples/MassAgeFehSamples_'+self.system+'.txt'

        self.phi_keys=[]
        for key,value in self.sampling_parameters.items():
            if value['dist']=='Samples':self.phi_keys.append(key)

        self.phi_mean=[]
        self.phi_sigma=[]
        for key in self.phi_keys:
            self.phi_mean=numpy.append(self.phi_mean,self.current_parameters[key])
            self.phi_sigma=numpy.append(self.phi_sigma,self.sampling_parameters[key]['step'])

        self.phi_sigma_matrix=numpy.zeros([len(self.phi_keys),len(self.phi_keys)])
        i=numpy.diag_indices(len(self.phi_keys))
        self.phi_sigma_matrix[i]=(1.0/self.phi_sigma**2)
        self.proposed=dict()

