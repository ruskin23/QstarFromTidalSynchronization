#!/usr/bin/env python3

from matplotlib import pyplot
from scipy import stats
import numpy


from dynesty import NestedSampler, plotting as dyplot

def prior_transform(unit_cube):
    """Insert fixd value between each input rv."""

    return unit_cube
    result = numpy.empty(unit_cube.size * 2 - 1)
    result[::2] = unit_cube
    result[1::2] = numpy.arange(unit_cube.size - 1)
    return result

def log_likelihood(parameters):
    """For now just gaussian mean=0.5, sigma=0.1 for the first param."""

    return numpy.log10(stats.norm.pdf(parameters[0], loc=0.5, scale=0.1))

def sample(nrandom):
    """Create and run the sampler return the resluts."""

    sampler = NestedSampler(log_likelihood,
                            prior_transform,
                            nlive=50000,
                            ndim=nrandom,#(2*nrandom - 1),
                            npdim=nrandom)
    sampler.run_nested(dlogz=0, maxcall=100000)
    return sampler.results

def corner_plots(results):
    """Create and display a corner plot."""

#    nparams = 2 * nrandom + 1
#    figure, axes = pyplot.subplots(nparams)
#    axes.reshape((nparams, nparams))

    dyplot.cornerplot(results, show_titles=True, quantiles=[0.025, 0.975])
    pyplot.show()

if __name__ == '__main__':
    corner_plots(sample(3))