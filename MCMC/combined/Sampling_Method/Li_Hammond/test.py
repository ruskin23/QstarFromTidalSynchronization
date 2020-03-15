from scipy import special
import numpy
import matplotlib.pyplot as plt
from utils import cummulative_distribution

sample_filename='../../mcmc_mass_age/samples/updated_samples/MassAgeFehSamples_10.txt'


mass_samples=[]
feh_samples=[]
multiplicity=[]

with open(sample_filename,'r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        mass_samples=numpy.append(mass_samples,float(x[0]))
        feh_samples=numpy.append(feh_samples,float(x[2]))
        multiplicity=numpy.append(multiplicity,float(x[3]))


CDF_mass=cummulative_distribution(mass_samples,multiplicity)()
CDF_feh=cummulative_distribution(feh_samples,multiplicity)()

plt.plot(*zip(*CDF_mass))


x=numpy.linspace(0,3,len(mass_samples))

f=(x-numpy.mean(mass_samples))/0.2
cdf=0.5*(1+special.erf(f/numpy.sqrt(2)))

print(CDF_mass[0:10])
print(cdf[0:10])

plt.plot(x,cdf)

plt.show()
