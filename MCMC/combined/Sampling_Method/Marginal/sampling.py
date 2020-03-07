import numpy
from covariance_matrix import Covariance
import random
import matplotlib.pyplot as plt
from utils import cummulative_distribution

class Sampling:

    def _multivariate_gaussian(self,
                               x_vector,
                               x_mean,
                               y_vector,
                               y_mean,
                               sigma_xy):

        arg=numpy.matmul(numpy.transpose(x_vector-x_mean),numpy.matmul(sigma_xy,(y_vector-y_mean)))
        if x_vector.all()!=y_vector.all():return numpy.exp(-2*arg)
        else:return numpy.exp(-arg)


    def stepping_function_normalization(self,
                                        phi_samples):


        return self.phi_summation(phi_samples)



    def joint_probability(self,
                          phi_N,
                          phi_mean,
                          theta_N,
                          theta_mean):


        phi_phi = self._multivariate_gaussian(phi_N,phi_mean,phi_N,phi_mean,self.SIGMAINV_phi_phi)
        theta_theta = self._multivariate_gaussian(theta_N,theta_mean,theta_N,theta_mean,self.SIGMAINV_theta_theta)
        theta_phi = self._multivariate_gaussian(theta_N,theta_mean,phi_N,phi_mean,self.SIGMAINV_theta_phi)

        return phi_phi*theta_theta*theta_phi

    def get_theta_samples(self):


        phi_N=self.phi_samples
        phi_mean=numpy.array([self.parameters[key] for key in self.phi_keys])
        theta_mean=numpy.array([self.parameters[key] for key in self.theta_keys])

        mean=theta_mean - numpy.matmul(self.SIGMA_theta_theta,numpy.matmul(self.SIGMAINV_theta_phi,numpy.transpose(phi_N-phi_mean)))
        cov=self.SIGMA_theta_theta

        return numpy.random.multivariate_normal(mean,cov)


    def phi_summation(self,phi_0):

        mean_vector=phi_0
        P=0
        with open(self.samples_file,'r') as f:
            next(f)
            for lines in f:
                x=lines.split()
                sample_vector=[]
                m_i=float(x[3])
                for keys in self.phi_keys:
                    sample_vector=numpy.append(sample_vector,float(x[self.phi_keys.index(keys)]))
                arg=numpy.matmul(numpy.transpose(sample_vector-mean_vector),numpy.matmul(self.SIGMAINV_phi_phi,(sample_vector-mean_vector)))
                P=P+m_i*numpy.exp(-arg)

        return P

    def get_phi_samples(self):


        mean_vector=numpy.array([self.parameters[key] for key in self.phi_keys])

        N=self.phi_summation(mean_vector)

        modified_multiplicity=[]
        Samples=[]
        with open(self.samples_file,'r') as f:
            next(f)
            for lines in f:
                x=lines.split()
                m_i=float(x[3])
                sample_vector=[]
                for keys in self.phi_keys:
                    sample_vector=numpy.append(sample_vector,float(x[self.phi_keys.index(keys)]))
                arg=numpy.matmul(numpy.transpose(sample_vector-mean_vector),numpy.matmul(self.SIGMAINV_phi_phi,(sample_vector-mean_vector)))
                modified_multiplicity=numpy.append(modified_multiplicity,(m_i*numpy.exp(-arg))/N)
                Samples.append(sample_vector)

        U=random.uniform(0,1)
        for i in range(len(modified_multiplicity)):
            if U<modified_multiplicity[i]:return Samples[i]
            else:U=U-modified_multiplicity[i]

    def __call__(self):

        self.phi_samples=self.get_phi_samples()
        self.theta_samples=self.get_theta_samples()


    def __init__(self,
                 system,
                 parameters
                 ):


        self.system=system
        self.parameters=parameters

        self.samples_file='../../../mcmc_mass_age/samples/updated_samples/MassAgeFehSamples_'+system+'.txt'

        self.phi_keys=['mass','age','feh']
        self.theta_keys=['Porb','e','Wdisk','logQ']

        C=Covariance(system)
        self.covariance_matrix=C.Calculate('Covariance')
        self.covariance_inverse=numpy.linalg.inv(self.covariance_matrix)

        n_theta=len(self.theta_keys)
        n_phi=len(self.phi_keys)

        self.SIGMA_theta_theta=self.covariance_matrix[0:n_theta,0:n_theta].copy()
        self.SIGMA_phi_phi=self.covariance_matrix[n_theta:n_theta+n_phi,n_theta:n_theta+n_phi].copy()
        self.SIGMA_theta_phi=self.covariance_matrix[0:n_theta,n_theta:n_theta+n_phi].copy()

        self.SIGMAINV_theta_theta=self.covariance_inverse[0:n_theta,0:n_theta].copy()
        self.SIGMAINV_phi_phi=self.covariance_inverse[n_theta:n_theta+n_phi,n_theta:n_theta+n_phi].copy()
        self.SIGMAINV_theta_phi=self.covariance_inverse[0:n_theta,n_theta:n_theta+n_phi].copy()




if __name__=='__main__':

    system='137'
    parameters=dict(mass=1.0,
                    age=4.0,
                    feh=-0.5,
                    Porb=8.0,
                    e=0.1,
                    Wdisk=4.1,
                    logQ=6)

    Sampling(system,parameters)()
