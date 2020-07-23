import argparse
import os
import sys
import warnings

import matplotlib.pyplot as plt

import numpy
from numpy import linalg
import scipy
from scipy import stats
import math
try:
    from scipy.special import logsumexp
except ImportError:
    from scipy.misc import logsumexp

import dynesty
from dynesty import plotting as dyplot
from dynesty import utils as dyfunc
from dynesty.dynamicsampler import weight_function
#from dynesty.utils import *
#from dynest.dynamicsampler import *


from pathos.pools import ProcessPool

import random
import copy

import pickle
import dill

from results import Results

def cmdline_args():

    parser = argparse.ArgumentParser()

    parser.add_argument('-l',
                        action='store',
                        dest='system',
                        help='select a system for mcmc'
                        )

    parser.add_argument('-n',
                        action='store',
                        dest='num_threads',
                        help='number of parallel processes'
                        )


    return parser.parse_args()

def _get_nsamps_samples_n(res):
    """ Helper function for calculating the number of samples
    Parameters
    ----------
    res : :class:`~dynesty.results.Results` instance
        The :class:`~dynesty.results.Results` instance taken from a previous
        nested sampling run.
    Returns
    -------
    nsamps: int
        The total number of samples
    samples_n: array
        Number of live points at a given iteration
    """
    try:
        # Check if the number of live points explicitly changes.
        samples_n = res.samples_n
        nsamps = len(samples_n)
    except:
        # If the number of live points is constant, compute `samples_n`.
        niter = res.niter
        nlive = res.nlive
        nsamps = len(res.logvol)
        if nsamps == niter:
            samples_n = numpy.ones(niter, dtype='int') * nlive
        elif nsamps == (niter + nlive):
            samples_n = numpy.append(numpy.ones(niter, dtype='int') * nlive,
                                  numpy.arange(1, nlive + 1)[::-1])
        else:
            raise ValueError("Final number of samples differs from number of "
                             "iterations and number of live points.")
    return nsamps, samples_n


def jitter_run(res, rstate=None, approx=False):
    """
    Probes **statistical uncertainties** on a nested sampling run by
    explicitly generating a *realization* of the prior volume associated
    with each sample (dead point). Companion function to :meth:`resample_run`
    and :meth:`simulate_run`.
    Parameters
    ----------
    res : :class:`~dynesty.results.Results` instance
        The :class:`~dynesty.results.Results` instance taken from a previous
        nested sampling run.
    rstate : `~numpy.random.RandomState`, optional
        `~numpy.random.RandomState` instance.
    approx : bool, optional
        Whether to approximate all sets of uniform order statistics by their
        associated marginals (from the Beta distribution). Default is `False`.
    Returns
    -------
    new_res : :class:`~dynesty.results.Results` instance
        A new :class:`~dynesty.results.Results` instance with corresponding
        weights based on our "jittered" prior volume realizations.
    """

    if rstate is None:
        rstate = numpy.random

    # Initialize evolution of live points over the course of the run.
    nsamps, samples_n = _get_nsamps_samples_n(res)
    logl = res.logl

    # Simulate the prior volume shrinkage associated with our set of "dead"
    # points. At each iteration, if the number of live points is constant or
    # increasing, our prior volume compresses by the maximum value of a set
    # of `K_i` uniformly distributed random numbers (i.e. as `Beta(K_i, 1)`).
    # If instead the number of live points is decreasing, that means we're
    # instead  sampling down a set of uniform random variables
    # (i.e. uniform order statistics).
    nlive_flag = numpy.ones(nsamps, dtype='bool')
    nlive_start, bounds = [], []

    if not approx:
        # Find all instances where the number of live points is either constant
        # or increasing.
        nlive_flag[1:] = numpy.diff(samples_n) >= 0

        # For all the portions that are decreasing, find out where they start,
        # where they end, and how many live points are present at that given
        # iteration.

        if numpy.any(~nlive_flag):
            i = 0
            while i < nsamps:
                if not nlive_flag[i]:
                    bound = []
                    bound.append(i-1)
                    nlive_start.append(samples_n[i-1])
                    while i < nsamps and not nlive_flag[i]:
                        i += 1
                    bound.append(i)
                    bounds.append(bound)
                i += 1

    # The maximum out of a set of `K_i` uniformly distributed random variables
    # has a marginal distribution of `Beta(K_i, 1)`.
    t_arr = numpy.zeros(nsamps)
    t_arr[nlive_flag] = rstate.beta(a=samples_n[nlive_flag], b=1)

    # If we instead are sampling the set of uniform order statistics,
    # we note that the jth largest value is marginally distributed as
    # `Beta(j, K_i-j+1)`. The full joint distribution is::
    #
    #     X_(j) / X_N = (Y_1 + ... + Y_j) / (Y_1 + ... + Y_{K+1})
    #
    # where X_(j) is the prior volume of the live point with the `j`-th
    # *best* likelihood (i.e. prior volume shrinks as likelihood increases)
    # and the `Y_i`'s are i.i.d. exponentially distributed random variables.
    nunif = len(nlive_start)
    for i in range(nunif):
        nstart = nlive_start[i]
        bound = bounds[i]
        sn = samples_n[bound[0]:bound[1]]
        y_arr = rstate.exponential(scale=1.0, size=nstart+1)
        ycsum = y_arr.cumsum()
        ycsum /= ycsum[-1]
        uorder = ycsum[numpy.append(nstart, sn-1)]
        rorder = uorder[1:] / uorder[:-1]
        t_arr[bound[0]:bound[1]] = rorder

    # These are the "compression factors" at each iteration. Let's turn
    # these into associated ln(volumes).
    logvol = numpy.log(t_arr).cumsum()

    # Compute weights using quadratic estimator.
    h = 0.
    logz = -1.e300
    loglstar = -1.e300
    logzvar = 0.
    logvols_pad = numpy.concatenate(([0.], logvol))
    logdvols = logsumexp(a=numpy.c_[logvols_pad[:-1], logvols_pad[1:]],
                         axis=1, b=numpy.c_[numpy.ones(nsamps), -numpy.ones(nsamps)])
    logdvols += math.log(0.5)
    dlvs = -numpy.diff(numpy.append(0., res.logvol))
    saved_logwt, saved_logz, saved_logzvar, saved_h = [], [], [], []
    for i in range(nsamps):
        loglstar_new = logl[i]
        logdvol, dlv = logdvols[i], dlvs[i]
        logwt = numpy.logaddexp(loglstar_new, loglstar) + logdvol
        logz_new = numpy.logaddexp(logz, logwt)
        lzterm = (math.exp(loglstar - logz_new) * loglstar +
                  math.exp(loglstar_new - logz_new) * loglstar_new)
        h_new = (math.exp(logdvol) * lzterm +
                 math.exp(logz - logz_new) * (h + logz) -
                 logz_new)
        dh = h_new - h
        h = h_new
        logz = logz_new
        logzvar += dh * dlv
        loglstar = loglstar_new
        saved_logwt.append(logwt)
        saved_logz.append(logz)
        saved_logzvar.append(logzvar)
        saved_h.append(h)

    # Copy results.
    new_res = Results([item for item in res.items()])

    # Overwrite items with our new estimates.
    new_res.logvol = numpy.array(logvol)
    new_res.logwt = numpy.array(saved_logwt)
    new_res.logz = numpy.array(saved_logz)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        new_res.logzerr = numpy.sqrt(numpy.array(saved_logzvar))
    new_res.h = numpy.array(saved_h)

    return new_res


def jitter_run(res, rstate=None, approx=False):
    """
    Probes **statistical uncertainties** on a nested sampling run by
    explicitly generating a *realization* of the prior volume associated
    with each sample (dead point). Companion function to :meth:`resample_run`
    and :meth:`simulate_run`.
    Parameters
    ----------
    res : :class:`~dynesty.results.Results` instance
        The :class:`~dynesty.results.Results` instance taken from a previous
        nested sampling run.
    rstate : `~numpy.random.RandomState`, optional
        `~numpy.random.RandomState` instance.
    approx : bool, optional
        Whether to approximate all sets of uniform order statistics by their
        associated marginals (from the Beta distribution). Default is `False`.
    Returns
    -------
    new_res : :class:`~dynesty.results.Results` instance
        A new :class:`~dynesty.results.Results` instance with corresponding
        weights based on our "jittered" prior volume realizations.
    """

    if rstate is None:
        rstate = numpy.random

    # Initialize evolution of live points over the course of the run.
    nsamps, samples_n = _get_nsamps_samples_n(res)
    logl = res.logl

    # Simulate the prior volume shrinkage associated with our set of "dead"
    # points. At each iteration, if the number of live points is constant or
    # increasing, our prior volume compresses by the maximum value of a set
    # of `K_i` uniformly distributed random numbers (i.e. as `Beta(K_i, 1)`).
    # If instead the number of live points is decreasing, that means we're
    # instead  sampling down a set of uniform random variables
    # (i.e. uniform order statistics).
    nlive_flag = numpy.ones(nsamps, dtype='bool')
    nlive_start, bounds = [], []

    if not approx:
        # Find all instances where the number of live points is either constant
        # or increasing.
        nlive_flag[1:] = numpy.diff(samples_n) >= 0

        # For all the portions that are decreasing, find out where they start,
        # where they end, and how many live points are present at that given
        # iteration.

        if numpy.any(~nlive_flag):
            i = 0
            while i < nsamps:
                if not nlive_flag[i]:
                    bound = []
                    bound.append(i-1)
                    nlive_start.append(samples_n[i-1])
                    while i < nsamps and not nlive_flag[i]:
                        i += 1
                    bound.append(i)
                    bounds.append(bound)
                i += 1

    # The maximum out of a set of `K_i` uniformly distributed random variables
    # has a marginal distribution of `Beta(K_i, 1)`.
    t_arr = numpy.zeros(nsamps)
    t_arr[nlive_flag] = rstate.beta(a=samples_n[nlive_flag], b=1)

    # If we instead are sampling the set of uniform order statistics,
    # we note that the jth largest value is marginally distributed as
    # `Beta(j, K_i-j+1)`. The full joint distribution is::
    #
    #     X_(j) / X_N = (Y_1 + ... + Y_j) / (Y_1 + ... + Y_{K+1})
    #
    # where X_(j) is the prior volume of the live point with the `j`-th
    # *best* likelihood (i.e. prior volume shrinks as likelihood increases)
    # and the `Y_i`'s are i.i.d. exponentially distributed random variables.
    nunif = len(nlive_start)
    for i in range(nunif):
        nstart = nlive_start[i]
        bound = bounds[i]
        sn = samples_n[bound[0]:bound[1]]
        y_arr = rstate.exponential(scale=1.0, size=nstart+1)
        ycsum = y_arr.cumsum()
        ycsum /= ycsum[-1]
        uorder = ycsum[numpy.append(nstart, sn-1)]
        rorder = uorder[1:] / uorder[:-1]
        t_arr[bound[0]:bound[1]] = rorder

    # These are the "compression factors" at each iteration. Let's turn
    # these into associated ln(volumes).
    logvol = numpy.log(t_arr).cumsum()

    # Compute weights using quadratic estimator.
    h = 0.
    logz = -1.e300
    loglstar = -1.e300
    logzvar = 0.
    logvols_pad = numpy.concatenate(([0.], logvol))
    logdvols = logsumexp(a=numpy.c_[logvols_pad[:-1], logvols_pad[1:]],
                         axis=1, b=numpy.c_[numpy.ones(nsamps), -numpy.ones(nsamps)])
    logdvols += math.log(0.5)
    dlvs = -numpy.diff(numpy.append(0., res.logvol))
    saved_logwt, saved_logz, saved_logzvar, saved_h = [], [], [], []
    for i in range(nsamps):
        loglstar_new = logl[i]
        logdvol, dlv = logdvols[i], dlvs[i]
        logwt = numpy.logaddexp(loglstar_new, loglstar) + logdvol
        logz_new = numpy.logaddexp(logz, logwt)
        lzterm = (math.exp(loglstar - logz_new) * loglstar +
                  math.exp(loglstar_new - logz_new) * loglstar_new)
        h_new = (math.exp(logdvol) * lzterm +
                 math.exp(logz - logz_new) * (h + logz) -
                 logz_new)
        dh = h_new - h
        h = h_new
        logz = logz_new
        logzvar += dh * dlv
        loglstar = loglstar_new
        saved_logwt.append(logwt)
        saved_logz.append(logz)
        saved_logzvar.append(logzvar)
        saved_h.append(h)

    # Copy results.
    new_res = Results([item for item in res.items()])

    # Overwrite items with our new estimates.
    new_res.logvol = numpy.array(logvol)
    new_res.logwt = numpy.array(saved_logwt)
    new_res.logz = numpy.array(saved_logz)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        new_res.logzerr = numpy.sqrt(numpy.array(saved_logzvar))
    new_res.h = numpy.array(saved_h)

    return new_res


def kld_error(res, error='simulate', rstate=None, return_new=False,
              approx=False):
    """
    Computes the `Kullback-Leibler (KL) divergence
    <https://en.wikipedia.org/wiki/Kullback-Leibler_divergence>`_ *from* the
    discrete probability distribution defined by `res` *to* the discrete
    probability distribution defined by a **realization** of `res`.
    Parameters
    ----------
    res : :class:`~dynesty.results.Results` instance
        :class:`~dynesty.results.Results` instance for the distribution we
        are computing the KL divergence *from*.
    error : {`'jitter'`, `'resample'`, `'simulate'`}, optional
        The error method employed, corresponding to :meth:`jitter_run`,
        :meth:`resample_run`, and :meth:`simulate_run`, respectively.
        Default is `'simulate'`.
    rstate : `~numpy.random.RandomState`, optional
        `~numpy.random.RandomState` instance.
    return_new : bool, optional
        Whether to return the realization of the run used to compute the
        KL divergence. Default is `False`.
    approx : bool, optional
        Whether to approximate all sets of uniform order statistics by their
        associated marginals (from the Beta distribution). Default is `False`.
    Returns
    -------
    kld : `~numpy.ndarray` with shape (nsamps,)
        The cumulative KL divergence defined *from* `res` *to* a
        random realization of `res`.
    new_res : :class:`~dynesty.results.Results` instance, optional
        The :class:`~dynesty.results.Results` instance corresponding to
        the random realization we computed the KL divergence *to*.
    """

    # Define our original importance weights.
    logp2 = res.logwt - res.logz[-1]

    # Compute a random realization of our run.
    if error == 'jitter':
        new_res = jitter_run(res, rstate=rstate, approx=approx)
    elif error == 'resample':
        new_res, samp_idx = resample_run(res, rstate=rstate, return_idx=True)
        logp2 = logp2[samp_idx]  # re-order our original results to match
    elif error == 'simulate':
        new_res, samp_idx = resample_run(res, rstate=rstate, return_idx=True)
        new_res = jitter_run(new_res)
        logp2 = logp2[samp_idx]  # re-order our original results to match
    else:
        raise ValueError("Input `'error'` option '{0}' is not valid."
                         .format(error))

    # Define our new importance weights.
    logp1 = new_res.logwt - new_res.logz[-1]

    # Compute the KL divergence.
    kld = numpy.cumsum(numpy.exp(logp1) * (logp1 - logp2))

    if return_new:
        return kld, new_res
    else:
        return kld


def _kld_error(args):
    """ Internal `pool.map`-friendly wrapper for :meth:`kld_error` used by
    :meth:`stopping_function`."""

    # Extract arguments.
    results, error, approx = args

    return kld_error(results, error, rstate=numpy.random, return_new=True,
                     approx=approx)

def stopping_function(results, args=None, rstate=None, M=None,
                      return_vals=False):
    """
    The default stopping function utilized by :class:`DynamicSampler`.
    Zipped parameters are passed to the function via :data:`args`.
    Assigns the run a stopping value based on a weighted average of the
    stopping values for the posterior and evidence::
        stop = pfrac * stop_post + (1.- pfrac) * stop_evid
    The evidence stopping value is based on the estimated evidence error
    (i.e. standard deviation) relative to a given threshold::
        stop_evid = evid_std / evid_thresh
    The posterior stopping value is based on the fractional error (i.e.
    standard deviation / mean) in the Kullback-Leibler (KL) divergence
    relative to a given threshold::
        stop_post = (kld_std / kld_mean) / post_thresh
    Estimates of the mean and standard deviation are computed using `n_mc`
    realizations of the input using a provided `'error'` keyword (either
    `'jitter'` or `'simulate'`, which call related functions :meth:`jitter_run`
    and :meth:`simulate_run` in :mod:`dynesty.utils`, respectively, or
    `'sim_approx'`, which boosts `'jitter'` by a factor of two).
    Returns the boolean `stop <= 1`. If `True`, the :class:`DynamicSampler`
    will stop adding new samples to our results.
    Parameters
    ----------
    results : :class:`Results` instance
        :class:`Results` instance.
    args : dictionary of keyword arguments, optional
        Arguments used to set the stopping values. Default values are
        `pfrac = 1.0`, `evid_thresh = 0.1`, `post_thresh = 0.02`,
        `n_mc = 128`, `error = 'sim_approx'`, and `approx = True`.
    rstate : `~numpy.random.RandomState`, optional
        `~numpy.random.RandomState` instance.
    M : `map` function, optional
        An alias to a `map`-like function. This allows users to pass
        functions from pools (e.g., `pool.map`) to compute realizations in
        parallel. By default the standard `map` function is used.
    return_vals : bool, optional
        Whether to return the stopping value (and its components). Default
        is `False`.
    Returns
    -------
    stop_flag : bool
        Boolean flag indicating whether we have passed the desired stopping
        criteria.
    stop_vals : tuple of shape (3,), optional
        The individual stopping values `(stop_post, stop_evid, stop)` used
        to determine the stopping criteria.
    """

    # Initialize values.
    if args is None:
        args = dict({})
    if rstate is None:
        rstate = numpy.random
    if M is None:
        M = map

    # Initialize hyperparameters.
    pfrac = args.get('pfrac', 1.0)
    if not 0. <= pfrac <= 1.:
        raise ValueError("The provided `pfrac` {0} is not between 0. and 1."
                         .format(pfrac))
    evid_thresh = args.get('evid_thresh', 0.1)
    if pfrac < 1. and evid_thresh < 0.:
        raise ValueError("The provided `evid_thresh` {0} is not non-negative "
                         "even though `1. - pfrac` is {1}."
                         .format(evid_thresh, 1. - pfrac))
    post_thresh = args.get('post_thresh', 0.02)
    if pfrac > 0. and post_thresh < 0.:
        raise ValueError("The provided `post_thresh` {0} is not non-negative "
                         "even though `pfrac` is {1}."
                         .format(post_thresh, pfrac))
    n_mc = args.get('n_mc', 128)
    if n_mc <= 1:
        raise ValueError("The number of realizations {0} must be greater "
                         "than 1.".format(n_mc))
    elif n_mc < 20:
        warnings.warn("Using a small number of realizations might result in "
                      "excessively noisy stopping value estimates.")
    error = args.get('error', 'sim_approx')
    if error not in {'jitter', 'simulate', 'sim_approx'}:
        raise ValueError("The chosen `'error'` option {0} is not valid."
                         .format(error))
    if error == 'sim_approx':
        error = 'jitter'
        boost = 2.
    else:
        boost = 1.
    approx = args.get('approx', True)

    # Compute realizations of ln(evidence) and the KL divergence.
    rlist = [results for i in range(n_mc)]
    error_list = [error for i in range(n_mc)]
    approx_list = [approx for i in range(n_mc)]
    args = zip(rlist, error_list, approx_list)
    outputs = list(M(_kld_error, args))
    kld_arr, lnz_arr = numpy.array([(kld[-1], res.logz[-1])
                                 for kld, res in outputs]).T

    # Evidence stopping value.
    lnz_std = numpy.std(lnz_arr)
    stop_evid = numpy.sqrt(boost) * lnz_std / evid_thresh

    # Posterior stopping value.
    kld_mean, kld_std = numpy.mean(kld_arr), numpy.std(kld_arr)
    stop_post = boost * (kld_std / kld_mean) / post_thresh

    # Effective stopping value.
    stop = pfrac * stop_post + (1. - pfrac) * stop_evid

    if return_vals:
        return stop <= 1., (stop_post, stop_evid, stop)
    else:
        return stop <= 1.


class DahiyaChutiya:

    def loglikelihood(self,x):

        # define likelihood constants
        ndim = 3
        C = numpy.identity(ndim)
        C[C==0] = 0.95
        Cinv = linalg.inv(C)
        lnorm = -0.5 * (numpy.log(2 * numpy.pi) * ndim + numpy.log(linalg.det(C)))

        y=numpy.array([x[4],x[5],x[6]])

        #for i in range(1000000):
        #    u=random.random()

        return -0.5 * numpy.dot(y, numpy.dot(Cinv, y)) + lnorm


    def prior_transform(self,u):

        x=numpy.array(u)

        for i,s in enumerate(self.sampling_parameters):

            if s[-1]=='Normal':
                mean=s[1]
                sigma=s[2]
                x[i]=scipy.stats.norm.ppf(u[i],loc=mean,scale=sigma)
            elif s[-1]=='Turncated_Normal':
                mean=s[1]
                sigma=s[2]
                low=(s[3]-mean)/sigma
                high=(s[4]-mean)/sigma
                x[i]=scipy.stats.truncnorm.ppf(u[i],low,high,loc=mean,scale=sigma)
            elif s[-1]=='Uniform':
                x[i]=(s[2]-s[1])*u[i] + s[1]

        return x


    def sample_initial(self,
                       dsampler,
                       print_func):

        live_points=None

        #Sample Initial Batch
        for results in dsampler.sample_initial(nlive=self.nlive_init,
                                               dlogz=self.dlogz_init,
                                               maxcall=self.maxcall_init,
                                               maxiter=self.maxiter_init,
                                               logl_max=self.logl_max_init,
                                               n_effective=self.n_effective_init,
                                               live_points=live_points):

            (worst, ustar, vstar, loglstar, logvol,
             logwt, logz, logzvar, h, nc, worst_it,
             boundidx, bounditer, eff, delta_logz) = results

            self.ncall+=nc
            self.niter+=1

            #print progress
            print_func(results,self.niter,self.ncall,nbatch=0,dlogz=self.dlogz_init,logl_max=self.logl_max_init)

        dsampler_dill=dsampler
        with open('sample_initial.dill','wb') as f:
            dill.dump(dsampler_dill,f)
            dill.dump(self.ncall,f)
            dill.dump(self.niter,f)

        return results




    def add_batches(self,
                    dsampler_dill,
                    print_func):


        for n in range(dsampler_dill.batch,self.maxbatch):

            print('\nBatch = {}'.format(dsampler_dill.batch))

            res=dsampler_dill.results
            mcall=min(self.maxcall-self.ncall,self.maxcall_batch)
            miter=min(self.maxiter-self.niter,self.maxiter_batch)
            neff=dsampler_dill.n_effective

            print('mcall = {}, miter = {}, neff = {}'.format(mcall,miter,neff))

            if mcall > 0  and miter > 0 and neff < numpy.inf:
                print('Evaluating Stopping Conditions')
                stop,stop_vals = stopping_function(res, self.stop_kwargs,
                                           rstate=dsampler_dill.rstate, M=dsampler_dill.M,
                                           return_vals=True)


                stop_post, stop_evid, stop_val = stop_vals
                print(stop)
            if mcall>0 and miter>0 and neff<numpy.inf and not stop:
                addBatch = dsampler_dill.add_batch(nlive=self.nlive_batch,
                                              wt_function=self.wt_function,
                                              wt_kwargs=self.wt_kwargs,
                                              maxiter=miter,
                                              maxcall=mcall,
                                              save_bounds=self.save_bounds,
                                              print_progress=self.print_progress,
                                              print_func=print_func,
                                              stop_val=stop_val)


                self.ncall,self.niter,logl_bounds,results=addBatch
                with open('sampler_at_batch_'+str(dsampler_dill.batch)+'.dill','wb') as f:
                    dill.dump(dsampler_dill,f)
                    dill.dump(self.ncall,f)
                    dill.dump(self.niter,f)
            elif logl_bounds[1]!=numpy.inf:
                break
            else:break





    def sampler(self):

        return dynesty.DynamicNestedSampler(self.loglikelihood, self.prior_transform,
                                       self.ndim,pool=self.pool,queue_size=self.queue_size)



    def __call__(self,
                 borof):


        if borof=='start':
            dsampler=self.sampler()
            print_func=None
            pbar,print_func=dsampler._get_print_func(print_func,self.print_progress)
            results=self.sample_initial(dsampler, print_func)

        elif borof=='continue':
            with open('sample_initial.dill','rb') as f:
                dsampler_dill=dill.load(f)
                self.ncall=dill.load(f)
                self.niter=dill.load(f)

            dsampler_dill.sampler.rstate=numpy.random
            pool=ProcessPool(nodes=4)
            dsampler_dill.sampler.pool=pool
            dsampler_dill.sampler.M=pool.map

            print_func=None
            pbar,print_func=dsampler_dill._get_print_func(print_func,self.print_progress)

            for n in range(dsampler_dill.batch,self.maxbatch):

                print('\nBatch = {}'.format(dsampler_dill.batch))

                res=dsampler_dill.results
                mcall=min(self.maxcall-self.ncall,self.maxcall_batch)
                miter=min(self.maxiter-self.niter,self.maxiter_batch)
                neff=dsampler_dill.n_effective

                print('mcall = {}, miter = {}, neff = {}'.format(mcall,miter,neff))

                if mcall > 0  and miter > 0 and neff < numpy.inf:
                    print('Evaluating Stopping Conditions')
                    stop,stop_vals = stopping_function(res, self.stop_kwargs,
                                               rstate=dsampler_dill.rstate, M=dsampler_dill.M,
                                               return_vals=True)


                    stop_post, stop_evid, stop_val = stop_vals
                    print(stop)
                if mcall>0 and miter>0 and neff<numpy.inf and not stop:
                    addBatch = dsampler_dill.add_batch(nlive=self.nlive_batch,
                                                  wt_function=self.wt_function,
                                                  wt_kwargs=self.wt_kwargs,
                                                  maxiter=miter,
                                                  maxcall=mcall,
                                                  save_bounds=self.save_bounds,
                                                  print_progress=self.print_progress,
                                                  print_func=print_func,
                                                  stop_val=stop_val)


                    self.ncall,self.niter,logl_bounds,results=addBatch
                    with open('sampler_at_batch_'+str(dsampler_dill.batch)+'.dill','wb') as f:
                        dill.dump(dsampler_dill,f)
                        dill.dump(self.ncall,f)
                        dill.dump(self.niter,f)
                elif logl_bounds[1]!=numpy.inf:
                    break
                else:break



        else: print('shamshya nahi')




    def __init__(self,
                 sampling_parameters,
                 pool,
                 queue_size):


        self.sampling_parameters=sampling_parameters
        self.pool=pool
        self.queue_size=queue_size

        # Initialize values.
        self.maxcall = sys.maxsize
        self.maxiter = sys.maxsize
        self.maxiter_batch = sys.maxsize
        self.maxcall_batch = sys.maxsize
        self.maxbatch = sys.maxsize
        self.maxiter_init = sys.maxsize
        self.maxcall_init = sys.maxsize
        self.wt_function = weight_function
        self.wt_kwargs = dict()

        self.stop_function = None
        self.stop_kwargs = dict()
        self.save_bounds=True
        self.M=pool.map

        #initialize call values
        self.ndim=len(self.sampling_parameters)
        self.nlive_init=500
        self.nlive_batch=500
        self.dlogz_init=0.01
        self.n_effective_init=numpy.inf
        self.logl_max_init=numpy.inf
        self.ncall=0
        self.it=1
        self.niter=self.it-1
        self.logl_bounds=(-numpy.inf,numpy.inf)
        self.maxcall_init = min(self.maxcall_init, self.maxcall)  # set max calls
        self.maxiter_init = min(self.maxiter_init, self.maxiter)  # set max iterations

        self.print_progress=True

args = cmdline_args()
system_number=args.system
num_threads=int(args.num_threads)
catalog_file='SpinlogQCatalog_el0.4.txt'

queue_size=num_threads
pool=ProcessPool(nodes=num_threads)

with open(catalog_file,'r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        at_system=x[0]
        if system_number==at_system:
            teff_value=float(x[2])
            teff_error=float(x[3])
            feh_value=float(x[4])
            feh_error=float(x[5])
            logg_value=float(x[10])
            logg_error=float(x[15])
            Porb_value=float(x[6])
            Porb_error=float(x[7])
            eccentricity_value=float(x[8])
            eccentricity_error=float(x[9])
            Pspin_value=float(x[12])
            Pspin_error=float(x[13])
            break



sampling_parameters = [('Porb',Porb_value,Porb_error,'Normal'),
                       ('feh',feh_value,feh_error,-1.014,0.537,'Truncated_normal'),
                       ('eccentricity',eccentricity_value,eccentricity_error,0.0,0.45,'Truncated_normal'),
                       ('Wdisk',2*scipy.pi/14,2*scipy.pi/1.4,'Uniform'),
                       ('logQ',5.0,12.0,'Uniform'),
                       ('primary_mass',0.5,1.2,'Uniform'),
                       ('age',1e-3,10.0,'Uniform')]



S=DahiyaChutiya(sampling_parameters,pool,queue_size)
S('continue')
