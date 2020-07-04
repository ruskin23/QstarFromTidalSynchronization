# import environment
from __future__ import division, print_function
from six.moves import range
import time, sys, os
import numpy as np
from numpy import linalg
import dynesty


# define likelihood constants
ndim = 3
C = np.identity(ndim)
C[C==0] = 0.95
Cinv = linalg.inv(C)
lnorm = -0.5 * (np.log(2 * np.pi) * ndim + np.log(linalg.det(C)))

# 3-D correlated multivariate normal log-likelihood
def loglikelihood(x):
    """Multivariate normal log-likelihood."""
    return -0.5 * np.dot(x, np.dot(Cinv, x)) + lnorm

# prior transform
def prior_transform(u):
    """Transforms our unit cube samples `u` to a flat prior between -10. and 10. in each variable."""
    return 10. * (2. * u - 1.)

# gradient of log-likelihood *with respect to u*
# -> d(lnl)/du = d(lnl)/dv * dv/du
# dv/du = 1. / prior(v)
def gradient(x):
    """Multivariate normal log-likelihood gradient."""
    dlnl_dv = -np.dot(Cinv, x)  # standard gradient
    jac = np.diag(np.full_like(x, 20.))  # Jacobian
    return np.dot(jac, dlnl_dv)  # transformed gradient


psampler = dynesty.NestedSampler(loglikelihood, prior_transform, ndim,
                                 nlive=1500, sample='hslice', gradient=gradient,
                                 first_update={'min_ncall': 5000, 'min_eff': 50.})

psampler.run_nested()
pres = psampler.results

