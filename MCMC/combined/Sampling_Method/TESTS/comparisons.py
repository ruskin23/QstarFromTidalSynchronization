import numpy
import matplotlib.pyploy as plt


class ComparisonPlots:


    def get_mcmc_distribution(self,
                              parameter):

        p=[]
        with open(self.mcmc_file,'r') as f:
            for lines in f:
                x=lines.split()
                if parameter=='age':p=numpy.append(p,float(x[6]))
                if parameter=='logQ':p=numpy.append(p,float(x[4]))

        self.sample_space=len(p)
        cdf=cummulative_distribution(p)
        return cdf()



    def __init__(system,
                 test_case):


        self.system=syetem
        self.sampling_method=sampling_method
        self.test_case=test_case

        self.sample_file='../../../mcmc_mass_age/samples/MassAgeFehSamples_'+self.system+'.txt'
        self.mcmc_file='../Uncorrelated_Phi/'+self.test_case+'/AcceptedParameters.txt'

        self.sample_space=0
        self.phi_mean_vector=
        self.model_width
        self.phi_sigminv_matrix=

