import numpy
from scipy import special
from scipy import integrate
import matplotlib.pyplot as plt
import itertools
import sys

class TestCorrelated:

    def _cummulative_distribution(self,
                                  x):

        a=[]
        v=0
        value=[]

        for tuples in x:
            a=numpy.append(a,tuples[0])
            v=v+tuples[1]
            value=numpy.append(value,v)

        return list(zip(a,value/max(value)))


    def __init__(self,
                 param_values):

        for key,value in param_values.items():
            setattr(self,key,value)

        self.mcmc_file='/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/combined/MCMC_TESTS/TEST5/AccetedParameters.txt'
        self.age_mcmc_samples=[]
        self.logQ_mcmc_samples=[]
        with open(self.mcmc_file,'r') as f:
            next(f)
            for lines in f:
                x=lines.split()
                self.age_mcmc_samples=numpy.append(self.age_mcmc_samples,float(x[6]))
                self.logQ_mcmc_samples=numpy.append(self.logQ_mcmc_samples,float(x[4]))


        self.samples_file='/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/mcmc_mass_age/samples/updated_samples/MassAgeFehSamples_54.txt'
        self.age_samples=[]
        self.mulitplicity=[]
        with open(self.samples_file,'r') as f:
            next(f)
            for lines in f:
                x=lines.split()
                self.age_samples=numpy.append(self.age_samples,float(x[1]))
                self.mulitplicity=numpy.append(self.mulitplicity,float(x[3]))
        self.logQ_samples=numpy.linspace(self.logQ_min,self.logQ_max,1000)



    def intergand(self,
                  t,
                  logQ):

        return  numpy.exp( - ((logQ-self.logQ_mean)/(self.sigma_logQ))**2 - (((t*logQ)-(self.t_mean*self.logQ_mean))/(self.sigma_t*self.sigma_logQ))**2 )


    def EvaluateProbability(self,
                            t,
                            m):

        I = lambda logQ: self.intergand(t,logQ)
        integral=integrate.quad(I,self.logQ_min,self.logQ_max)[0]
        return m*numpy.exp(-((t-self.t_mean)**2)/(self.sigma_t**2))*integral


    def CompareAge(self):

        N=0
        modified_multiplicity=[]
        for t,m in zip(self.age_samples,self.mulitplicity):
            P=self.EvaluateProbability(t,m)
            modified_multiplicity=numpy.append(modified_multiplicity,P)
            N=N+P

        modified_multiplicity=modified_multiplicity/N
        age_tuple=list(zip(self.age_samples,modified_multiplicity))
        age_tuple=sorted(age_tuple, key=lambda tup: tup[0])
        age_cummulative=self._cummulative_distribution(age_tuple)
        plt.scatter(*zip(*age_cummulative),color='r',label='analytical')

        age_mcmc_tuple=sorted([(x, len(list(y))) for x, y in itertools.groupby(self.age_mcmc_samples)], key=lambda tup: tup[0])
        age_mcmc_cummulative=self._cummulative_distribution(age_mcmc_tuple)
        plt.plot(*zip(*age_mcmc_cummulative),label='MCMC')

        plt.legend()
        plt.show()


    def CompareLogQ(self):

        Sum=numpy.zeros(len(self.logQ_samples))
        for t,m in zip(self.age_samples,self.mulitplicity):
            arg = ((t-self.t_mean)/self.sigma_t)**2 + ((t*self.logQ_samples - self.t_mean*self.logQ_mean)/(self.sigma_t*self.sigma_logQ))**2
            Sum = Sum + m*numpy.exp(-arg)

        logQ_dist=numpy.exp(-((self.logQ_samples-self.logQ_mean)/self.sigma_logQ)**2)*Sum
        logQ_tuple=list(zip(self.logQ_samples,logQ_dist))
        logQ_cummulative=self._cummulative_distribution(logQ_tuple)
        plt.scatter(*zip(*logQ_cummulative),color='r',label='analytical')


        logQ_mcmc_tuple=[(x, len(list(y))) for x, y in itertools.groupby(self.logQ_mcmc_samples)]
        logQ_mcmc_tuple=sorted(logQ_mcmc_tuple, key=lambda tup: tup[0])
        logQ_mcmc_cummulative=self._cummulative_distribution(logQ_mcmc_tuple)
        plt.plot(*zip(*logQ_mcmc_cummulative),label='MCMC')

        plt.legend()
        plt.show()



    def __call__(self,
                 param):

        if param=='logQ':self.CompareLogQ()
        if param=='age':self.CompareAge()

if __name__=='__main__':


    param_values=dict(logQ_mean=7.125744153362705,
                      t_mean=4.049835253026934,

                      sigma_logQ=0.2,
                      sigma_t=0.8,

                      logQ_min=5.0,
                      logQ_max=12.0
                      )


    Compare=TestCorrelated(param_values)
    Compare(sys.argv[1])

