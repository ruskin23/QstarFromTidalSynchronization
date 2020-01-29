import scipy
from scipy.stats import norm

import numpy
import matplolib.pyplot as plt



parameter=dict(value=1.0,
        sigma=0.5)

a=[]
for i in range(10000):
    a.append(scipy.stats.norm.rvs(loc=parameters['value'],scale=parameters['sigma']))

x=numpy.linspace(-0.01,0.01,len(a))

plt.plot(x,a)
pltshow()

