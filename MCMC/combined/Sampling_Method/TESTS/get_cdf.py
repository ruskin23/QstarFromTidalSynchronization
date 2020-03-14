import numpy
import sys
sys.path.append('../../')
from utils import cummulative_distribution
from test_models import TestModels

class GetCDF:

    def get_mcmc_distribution(self,
                              test_case,
                              parameter):

        mcmc_file='../Adaptive/'+test_case+'/AcceptedParameters.txt'
        p=[]
        with open(mcmc_file,'r') as f:
            for lines in f:
                x=lines.split()
                if parameter=='age':p=numpy.append(p,float(x[6]))
                if parameter=='logQ':p=numpy.append(p,float(x[4]))

        cdf=cummulative_distribution(p)

        return cdf()

    def get_model_distribution(self,
                               test_case,
                               parameter=None,
                               phi=None,
                               theta=None):



        if test_case=='prior':return self.Model.prior_model(phi,parameter)
        if test_case=='gp':return self.Model.gaussian_phi_model(phi)
        if test_case=='gpt':return self.Model.gaussian_phi_theta_model(parameter,phi)
        if test_case=='gptc':return self.Model.correlated_phi_theta(parameter,phi)
        if test_case=='gtt':return self.Model.gaussian_theta_theta_model(theta,parameter)
        if test_case=='gttc':return self.Model.correlated_theta_theta(theta,parameter)


    def __init__(self,
                 system,
                 sampling_parameters,
                 model_width,
                 rho):


        self.system=system
        self.sampling_parameters=sampling_parameters

        self.sample_file='../../../mcmc_mass_age/samples/updated_samples/MassAgeFehSamples_'+self.system+'.txt'

        self.Model=TestModels(self.sample_file,
                              sampling_parameters,
                              model_width,
                              rho)
