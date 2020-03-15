import matplotlib.pyplot as plt
import numpy
import sys
import itertools

class PercentileAge:

    def _cummulative_distribution(self,
                                  x):
        a=[]
        v=0.0
        values=[]

        for tuples in x:
            a=numpy.append(a,tuples[0])
            v=v+tuples[1]
            values=numpy.append(values,v)
        return list(zip(a,values/max(values)))

    def __init__(self,
                 system,
                 percentile):
        self.system=system
        self.percentile=percentile/100
        self.samplefile='/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/mcmc_mass_age/samples/updated_samples/MassAgeFehSamples_'+system+'.txt'

    def __call__(self):

        ages=[]
        multiplicity=[]
        with open(self.samplefile,'r') as f:
            next(f)
            for lines in f:
                x=lines.split()
                ages=numpy.append(ages,float(x[1]))
                multiplicity=numpy.append(multiplicity,float(x[3]))

        ages_tuples=zip(ages,multiplicity)
        ages_tuples=sorted(ages_tuples, key=lambda tup: tup[0])
        ages_cummulative=self._cummulative_distribution(ages_tuples)
        for age in ages_cummulative:
            if abs(self.percentile-age[1])<0.01:return age[0]
