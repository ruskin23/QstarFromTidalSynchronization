import numpy
import scipy
from scipy.stats import truncnorm
from scipy.stats import norm
import matplotlib.pyplot as plt


lower, higher = 0.1,2
mu,sigma = 1,0.5

X=scipy.stats.truncnorm(lower,higher,loc=mu,scale=sigma)
N=scipy.stats.norm(loc=mu,scale=sigma)


fig, ax = plt.subplots(2)

ax[0].hist(X.rvs(10000), normed=True)
ax[1].hist(N.rvs(10000), normed=True)
plt.show()

