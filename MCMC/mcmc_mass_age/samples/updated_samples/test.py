import numpy
import sys
import matplotlib.pyplot as plt
from scipy.stats import norm

class CDF:

    def sort_tuples(self):

        self.dist=sorted(self.dist, key=lambda tup: tup[0])

    def cummulative_distribution(self):

        self.sort_tuples()

        a=[]
        v=0
        value=[]

        for tuples in self.dist:
            a=numpy.append(a,tuples[0])
            v=v+tuples[1]
            value=numpy.append(value,v)

        return list(zip(a,value/max(value)))

    def __init__(self,
                 x):
        self.dist=x

system=sys.argv[1]
age_mean=7.249465124531277
age_sigma=0.5

age=[]
multiplicity=[]
modified_multiplicity=[]
no_multiplicity=[]
N=0
N1=0
with open('MassAgeFehSamples_'+system+'.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        age=numpy.append(age,float(x[1]))

        multiplicity=numpy.append(multiplicity,float(x[3]))
        modified_multiplicity=numpy.append(modified_multiplicity,float(x[3])*numpy.exp(-0.5*((float(x[1])-age_mean)/(age_sigma))**2))
        no_multiplicity=numpy.append(no_multiplicity,numpy.exp(-0.5*((float(x[1])-age_mean)/(age_sigma))**2))

        N1=N1+numpy.exp(-0.5*((float(x[1])-age_mean)/(age_sigma))**2)
        N=N+float(x[3])*numpy.exp(-0.5*((float(x[1])-age_mean)/(age_sigma))**2)

multiplicity=multiplicity/max(multiplicity)
modified_multiplicity=modified_multiplicity/N
#modified_multiplicity=modified_multiplicity/max(modified_multiplicity)
no_multiplicity=no_multiplicity/N1
#no_multiplicity=no_multiplicity/max(no_multiplicity)
cdf_m=CDF(list(zip(age,multiplicity))).cummulative_distribution()

cdf_modiified_multiplicity=CDF(list(zip(age,modified_multiplicity))).cummulative_distribution()
cdf_no_mutiplicity=CDF(list(zip(age,no_multiplicity))).cummulative_distribution()

#inverse_cdf_modiified_multiplicity=norm.ppf(cdf_modiified_multiplicity)
#inverse_cdf_no_mutiplicity=norm.ppf(cdf_no_mutiplicity)

#plt.scatter(age,multiplicity)
plt.plot(*zip(*cdf_m))
#plt.plot(*zip(*inverse_cdf_no_mutiplicity))
plt.show()
