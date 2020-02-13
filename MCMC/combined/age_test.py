import numpy
from numpy import inf
import matplotlib.pyplot as plt
import itertools

age_mean=4.049835253026934
sigma_theory=0.8
sigma_step=1.0


def norm(x,
         loc=0.0,
         sigma=1.0):

    y=x-loc
    y=y/sigma
    return numpy.exp(-(y**2)/2)


samples='/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/mcmc_mass_age/samples/MassAgeFehSamples_54.txt'
reduced_samples='/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/mcmc_mass_age/samples/updated_samples/MassAgeFehSamples_54.txt'
mcmc_file='/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/combined/kartof/MCMC_54/TEST2/AccetedParameters.txt'

age=[]
multiplicity=[]



with open(reduced_samples,'r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        age_i=float(x[1])
        m_i=float(x[3])
        age=numpy.append(age,age_i)
        multiplicity=numpy.append(multiplicity,m_i*norm(age_i-age_mean,loc=0,sigma=sigma_theory))


age_tuple=list(zip(age,multiplicity))
age_tuple=sorted(age_tuple, key=lambda tup: tup[0])

a=[]
value=[]
v=0
for p in age_tuple:
    a=numpy.append(a,p[0])
    v=v+p[1]
    value=numpy.append(value,v)

plt.scatter(a,value/max(value),color='r')
plt.axvline(x=age_mean)
#plt.show()

mcmc_age=[]
with open(mcmc_file,'r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        mcmc_age=numpy.append(mcmc_age,float(x[6]))

mcmc_age=[(x, len(list(y))) for x, y in itertools.groupby(mcmc_age)]
mcmc_age=sorted(mcmc_age, key=lambda tup: tup[0])

value=[]
multiplicity=[]
m=0
for a in mcmc_age:
    value=numpy.append(value,a[0])
    m=m+a[1]
    multiplicity=numpy.append(multiplicity,m)

multiplicity=multiplicity/max(multiplicity)
plt.plot(value,multiplicity)
plt.show()


