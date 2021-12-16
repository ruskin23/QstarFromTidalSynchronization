import numpy
from scipy import special

def erf_fun(x):
    return 0.5*(1+special.erf(x/numpy.sqrt(2)))
