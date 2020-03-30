import numpy
from covariance_matrix import Covariance
import random
import matplotlib.pyplot as plt
import sys
sys.path.append('../../')
from utils import cummulative_distribution
import scipy
from scipy import integrate
from scipy import special

class AdaptiveSampling:

    def _multivariate_gaussian(self,
                               x_vector,
                               x_mean,
                               y_vector,
                               y_mean,
                               sigma_xy):

        arg=numpy.matmul(numpy.transpose(x_vector-x_mean),numpy.matmul(sigma_xy,(y_vector-y_mean)))
        if x_vector.all()!=y_vector.all():return numpy.exp(-2*arg)
        else:return numpy.exp(-arg)

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
                arg=numpy.matmul(numpy.transpose(sample_vector-mean_vector),numpy.matmul(numpy.linalg.inv(self.SIGMA_phi_phi),(sample_vector-mean_vector)))
                P=P+m_i*numpy.exp(-0.5*arg)

        return P


    def integration(self,theta_mean_vector):
        x0=theta_mean_vector[1]
        y0=theta_mean_vector[2]
        z0=theta_mean_vector[3]
        return (-0.029099*special.erf(6.69861 - 0.558217*z0)*special.erf(70.0808 - 70.0238*x0
                                                                      + 0.0348903*y0 -
                                                                      0.00602739*z0)*special.erf(0.23186
                                                                                                 -
                                                                                                 0.611938*y0
                                                                                                 +
                                                                                                 0.00311609*z0) +
             0.029099*special.erf(6.69861 - 0.558217 *z0)*special.erf(0.0569769 - 70.0238
                                                                      *x0 + 0.0348903 *y0 -
                                                                      0.00602739
                                                                      *z0)*special.erf(0.23186
                                                                                       -
                                                                                       0.611938
                                                                                       *y0 +
                                                                                       0.00311609
                                                                                       *z0) +
             0.029099*special.erf(2.79109 - 0.558217 *z0)*special.erf(70.0386 - 70.0238*x0 +
                                                                      0.0348903*y0 -
                                                                      0.00602739*z0)*special.erf(0.253672
                                                                                                 -
                                                                                                 0.611938*y0
                                                                                                 +
                                                                                                 0.00311609*z0) -
             0.029099*special.erf(2.79109 - 0.558217*z0)*special.erf(0.0147852 - 70.0238*x0
                                                                     + 0.0348903*y0 -
                                                                     0.00602739*z0)*special.erf(0.253672
                                                                                                -
                                                                                                0.611938*y0
                                                                                                +
                                                                                                0.00311609*z0) +
             0.029099*special.erf(6.69861 - 0.558217*z0)*special.erf(69.9398 - 70.0238*x0 +
                                                                     0.0348903*y0 -
                                                                     0.00602739*z0)*special.erf(2.70409
                                                                                                -
                                                                                                0.611938*y0
                                                                                                +
                                                                                                0.00311609*z0) -
             0.029099*special.erf(2.79109 - 0.558217*z0)*special.erf(69.8976 - 70.0238*x0 +
                                                                     0.0348903*y0 -
                                                                     0.00602739*z0)*special.erf(2.7259
                                                                                                -
                                                                                                0.611938*y0
                                                                                                +
                                                                                                0.00311609*z0) +
             0.029099*special.erf(6.69861 - 0.558217*z0)*special.erf(2.70409 - 0.611938*y0 +
                                                                     0.00311609*z0)*special.erf(0.0839797
                                                                                                +
                                                                                                70.0238
                                                                                                *x0
                                                                                                -
                                                                                                0.0348903*y0
                                                                                                +
                                                                                                0.00602739*z0) -
             0.029099*special.erf(2.79109 - 0.558217*z0)*special.erf(2.7259 - 0.611938*y0 +
                                                                     0.00311609*z0)*special.erf(0.126171
                                                                                                +
                                                                                                70.0238
                                                                                                *x0
                                                                                                -
                                                                                                0.0348903*y0
                                                                                                +
                                                                                                0.00602739*z0))


    def stepping_function_normalization(self,
                                        phi_samples,
                                        theta_samples):

        phi_mean_vector=phi_samples
        theta_mean_vector=theta_samples

        P=0
        with open(self.samples_file,'r') as f:
            next(f)
            for lines in f:

                x=lines.split()
                phi_sample_vector=[]
                m_i=float(x[3])

                for keys in self.phi_keys:
                    phi_sample_vector=numpy.append(phi_sample_vector,float(x[self.phi_keys.index(keys)]))

                arg1=numpy.matmul(numpy.transpose(phi_sample_vector-phi_mean_vector),numpy.matmul(numpy.linalg.inv(self.SIGMA_phi_phi),(phi_sample_vector-phi_mean_vector)))

                V=numpy.matmul(self.V_theta_phi.T,numpy.matmul(numpy.linalg.inv(self.V_theta_theta),self.V_theta_phi))
                arg2=numpy.matmul(numpy.transpose(phi_sample_vector-phi_mean_vector),numpy.matmul(V,(phi_sample_vector-phi_mean_vector)))

                theta_mean_prime_vector=theta_mean_vector-numpy.matmul(numpy.linalg.inv(self.V_theta_theta),numpy.matmul(self.V_theta_phi,(phi_sample_vector-phi_mean_vector)))

                P=P+m_i*numpy.exp(-0.5*arg1)*numpy.exp(0.5*arg2)*self.integration(theta_mean_prime_vector)

        return P

    def get_theta_samples(self):


        phi_N=self.phi_samples
        print('Phi Samples = ',phi_N)
        phi_mean=numpy.array([self.parameters[key] for key in self.phi_keys])
        print('phi mean = ',phi_mean)
        theta_mean=numpy.array([self.parameters[key] for key in self.theta_keys])

        mean=theta_mean - numpy.matmul(numpy.linalg.inv(self.V_theta_theta),numpy.matmul(self.V_theta_phi,phi_N-phi_mean))
        cov=numpy.linalg.inv(self.V_theta_theta)

        print('MEAN = ',mean)

        samples=numpy.random.multivariate_normal(mean,cov)

        for i,key in enumerate(self.theta_keys):
            self.proposed[key]=samples[i]



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
                arg=numpy.matmul(numpy.transpose(sample_vector-mean_vector),numpy.matmul(numpy.linalg.inv(self.SIGMA_phi_phi),(sample_vector-mean_vector)))
                modified_multiplicity=numpy.append(modified_multiplicity,(m_i*numpy.exp(-0.5*arg))/N)
                Samples.append(sample_vector)

        U=random.uniform(0,1)
        for mm,s in zip(modified_multiplicity,Samples):
            if U<mm:return s
            else:U=U-mm

    def __call__(self):


        self.phi_samples=self.get_phi_samples()
        self.get_theta_samples()
        for k,keys in enumerate(self.phi_keys):
            self.proposed[keys]=self.phi_samples[k]
        return self.proposed

    def __init__(self,
                 system,
                 sampling_parameters,
                 parameters,
                 covariance_matrix
                 ):


        self.system=system
        self.parameters=parameters
        self.sampling_parameters=sampling_parameters

        self.samples_file='/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/mcmc_mass_age/samples/updated_samples/MassAgeFehSamples_'+system+'.txt'

        self.phi_keys=['primary_mass','age','feh']
        self.theta_keys=['Porb','eccentricity','Wdisk','logQ']

        self.covariance_matrix=covariance_matrix
        self.covariance_inverse=numpy.linalg.inv(self.covariance_matrix)

        n_theta=len(self.theta_keys)
        n_phi=len(self.phi_keys)

        self.SIGMA_theta_theta=self.covariance_matrix[0:n_theta,0:n_theta].copy()
        self.SIGMA_phi_phi=self.covariance_matrix[n_theta:n_theta+n_phi,n_theta:n_theta+n_phi].copy()
        self.SIGMA_theta_phi=self.covariance_matrix[0:n_theta,n_theta:n_theta+n_phi].copy()

        self.V_theta_theta=self.covariance_inverse[0:n_theta,0:n_theta].copy()
        self.V_phi_phi=self.covariance_inverse[n_theta:n_theta+n_phi,n_theta:n_theta+n_phi].copy()
        self.V_theta_phi=self.covariance_inverse[0:n_theta,n_theta:n_theta+n_phi].copy()

        self.proposed=dict()
