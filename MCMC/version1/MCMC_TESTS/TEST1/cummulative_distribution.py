import numpy
import matplotlib.pyplot as plt

data=[]

with open('AccetedParameters.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        data.append(float(x[6]))


values,base=numpy.histogram(data,bins=1000)
cummulative=numpy.cumsum(values)
cummulative=cummulative/max(cummulative)
plt.plot(base[:-1],cummulative,c='blue')

data=[]

with open('/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/mcmc_mass_age/samples/MassAgeFehSamples_54.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        data.append(float(x[2]))


values,base=numpy.histogram(data,bins=1000)
cummulative=numpy.cumsum(values)
cummulative=cummulative/max(cummulative)
plt.plot(base[:-1],cummulative,c='green')


plt.show()






