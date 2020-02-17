#correlation between mass and age
import sys
import numpy
from numpy import inf
import matplotlib.pyplot as plt
import itertools


def norm(x,
         loc=0.0,
         sigma=1.0):

    y=x-loc
    y=y/sigma
    return numpy.exp(-(y**2)/2)


param=sys.argv[1]

parameters=dict(mass=dict(mean=1.024563030621552,
                          sigma_theory=0.8,
                          sigma_step=1.0,),
                age=dict(mean=4.049835253026934,
                          sigma_theory=0.8,
                          sigma_step=1.0,),
                feh=dict(mean=-0.0508939712472557,
                          sigma_theory=0.2,
                          sigma_step=0.3,),
                logQ=dict(mean=7.125744153362705,
                          sigma_theory=0.2,
                          sigma_step=0.5)

                )

param_keys=['mass','age','feh']

samples='/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/mcmc_mass_age/samples/MassAgeFehSamples_54.txt'
reduced_samples='/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/mcmc_mass_age/samples/updated_samples/MassAgeFehSamples_54.txt'
mcmc_file='/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/combined/MCMC_TESTS/TEST5/AccetedParameters.txt'

multiplicity=[]
param_values=[]

with open(reduced_samples,'r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        mass_i=float(x[0])
        age_i=float(x[1])
        feh_i=float(x[2])
        m_i=float(x[3])

        if param=='age':
            param_values=numpy.append(param_values,age_i)
        if param=='mass':
            param_values=numpy.append(param_values,mass_i)
        if param=='feh':
            param_values=numpy.append(param_values,feh_i)

        factor1=norm(mass_i-parameters['mass']['mean'],loc=0,sigma=parameters['mass']['sigma_theory'])
        factor2=norm(age_i-parameters['age']['mean'],loc=0,sigma=parameters['age']['sigma_theory'])
        factor3=norm(mass_i*age_i-parameters['mass']['mean']*parameters['age']['mean'],loc=0,sigma=parameters['mass']['sigma_theory']*parameters['age']['sigma_theory'])
        factor=factor1*factor2*factor3

        multiplicity=numpy.append(multiplicity,m_i*factor)


param_tuple=list(zip(param_values,multiplicity))
param_tuple=sorted(param_tuple, key=lambda tup: tup[0])

a=[]
value=[]
v=0
for p in param_tuple:
    a=numpy.append(a,p[0])
    v=v+p[1]
    value=numpy.append(value,v)

plt.scatter(a,value/max(value),color='r')
plt.axvline(x=parameters[param]['mean'])



mcmc_param_values=[]
with open(mcmc_file,'r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        mcmc_param_values=numpy.append(mcmc_param_values,float(x[param_keys.index(param)+5]))

mcmc_param_values=[(x, len(list(y))) for x, y in itertools.groupby(mcmc_param_values)]
mcmc_param_values=sorted(mcmc_param_values, key=lambda tup: tup[0])

value=[]
multiplicity=[]
m=0
for a in mcmc_param_values:
    value=numpy.append(value,a[0])
    m=m+a[1]
    multiplicity=numpy.append(multiplicity,m)

multiplicity=multiplicity/max(multiplicity)
plt.plot(value,multiplicity)
plt.show()


