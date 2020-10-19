import numpy as np
import emcee
import time
def log_prob(x,mu,cov):
    diff=x-mu
    time.sleep(1)
    return -0.5 * np.dot(diff, np.linalg.solve(cov, diff))


ndim = 5

np.random.seed(42)
means = np.random.rand(ndim)

cov = 0.5 - np.random.rand(ndim ** 2).reshape((ndim, ndim))
cov = np.triu(cov)
cov += cov.T - np.diag(cov.diagonal())#pylint: disable=no-member
cov = np.dot(cov, cov)


nwalkers=20
p0=np.random.rand(nwalkers,ndim)

sampler=emcee.EnsembleSampler(nwalkers,ndim,log_prob,args=[means,cov])
state = sampler.run_mcmc(p0, 100,progress=True)