# -*- coding: utf-8 -*-

"""Modified HDF5 EMCEE backend EMCEE that buffers steps in case saving fails."""

from __future__ import division, print_function

__all__ = ["HDFBackend", "TempHDFBackend"]

import os
import os.path
from tempfile import NamedTemporaryFile
import pickle
import logging
from time import sleep

import numpy as np

from emcee import __version__
from emcee.backends import Backend


try:
    import h5py
except ImportError:
    h5py = None


class HDFBackend(Backend):
    """A backend that stores the chain in an HDF5 file using h5py
    .. note:: You must install `h5py <http://www.h5py.org/>`_ to use this
        backend.
    Args:
        filename (str): The name of the HDF5 file where the chain will be
            saved.
        name (str; optional): The name of the group where the chain will
            be saved.
        read_only (bool; optional): If ``True``, the backend will throw a
            ``RuntimeError`` if the file is opened with write access.
    """

    #Inherited from EMCEE
    #pylint: disable=invalid-name
    def _save_step_to_file(self, state, accepted, f):
        """Save a single step to already appropriately open HDF5 file."""

        g = f[self.name]
        iteration = g.attrs["iteration"]

        g["chain"][iteration, :, :] = state.coords
        g["log_prob"][iteration, :] = state.log_prob
        if state.blobs is not None:
            g["blobs"][iteration, :] = state.blobs
        g["accepted"][:] += accepted

        for i, v in enumerate(state.random_state):
            g.attrs["random_state_{0}".format(i)] = v

        g.attrs["iteration"] = iteration + 1
    #pylint: enable=invalid-name


    def _flush_unsaved_steps(self):
        """Save any unsaved steps and initialize a fresh unsaved steps file."""

        pending_steps = os.path.exists(self.unsaved_steps_fname)
        saved_iterations = 0
        if os.path.exists(self.filename):
            with self.open('r+' if pending_steps else 'r') as progress_file:
                if self.name not in progress_file:
                    return
                saved_iterations = progress_file[self.name].attrs['iteration']
                if pending_steps:
                    with open(self.unsaved_steps_fname, 'rb') as \
                            unsaved_steps_file:
                        unsaved_iteration = pickle.load(unsaved_steps_file)
                        assert unsaved_iteration <= saved_iterations
                        try:
                            while unsaved_iteration < saved_iterations:
                                pickle.load(unsaved_steps_file)
                                unsaved_iteration += 1
                            while True:
                                step = pickle.load(unsaved_steps_file)
                                self._save_step_to_file(*step, progress_file)
                                saved_iterations += 1
                        except EOFError:
                            logging.getLogger(__name__).info(
                                'Successfully added %d steps to %s',
                                saved_iterations - unsaved_iteration,
                                self.filename
                            )
        elif pending_steps:
            with open(self.unsaved_steps_fname, 'rb') as unsaved_steps_file:
                assert pickle.load(unsaved_steps_file) == 0

        with open(self.unsaved_steps_fname, 'wb') as unsaved_steps_file:
            pickle.dump(saved_iterations, unsaved_steps_file)


    #Inherited from EMCEE
    #pylint: disable=super-init-not-called
    def __init__(self, filename, name="mcmc", read_only=False, dtype=None):
        if h5py is None:
            raise ImportError("you must install 'h5py' to use the HDFBackend")
        self.filename = filename
        self.name = name
        self.read_only = read_only
        if dtype is None:
            self.dtype_set = False
            self.dtype = np.float64
        else:
            self.dtype_set = True
            self.dtype = dtype

        self._has_blobs = None
        self.unsaved_steps_fname = (os.path.splitext(filename)[0]
                                    +
                                    '.unsaved_steps')
        self._flush_unsaved_steps()
    #pylint: enable=super-init-not-called

    #Inherited from EMCEE
    #pylint: disable=missing-function-docstring
    #pylint: disable=invalid-name
    @property
    def initialized(self):
        if not os.path.exists(self.filename):
            return False
        try:
            with self.open() as f:
                return self.name in f
        except (OSError, IOError):
            return False

    def open(self, mode="r"):
        if self.read_only and mode != "r":
            raise RuntimeError(
                "The backend has been loaded in read-only "
                "mode. Set `read_only = False` to make "
                "changes."
            )
        f = h5py.File(self.filename, mode)
        if not self.dtype_set and self.name in f:
            g = f[self.name]
            if "chain" in g:
                self.dtype = g["chain"].dtype
                self.dtype_set = True
            self._has_blobs = g.attrs["has_blobs"]
        return f

    def reset(self, nwalkers, ndim):
        """Clear the state of the chain and empty the backend
        Args:
            nwakers (int): The size of the ensemble
            ndim (int): The number of dimensions
        """
        with self.open("a") as f:
            if self.name in f:
                del f[self.name]

            self._has_blobs = False

            g = f.create_group(self.name)
            g.attrs["version"] = __version__
            g.attrs["nwalkers"] = nwalkers
            g.attrs["ndim"] = ndim
            g.attrs["has_blobs"] = False
            g.attrs["iteration"] = 0
            g.create_dataset("accepted", data=np.zeros(nwalkers))
            g.create_dataset(
                "chain",
                (0, nwalkers, ndim),
                maxshape=(None, nwalkers, ndim),
                dtype=self.dtype,
            )
            g.create_dataset(
                "log_prob",
                (0, nwalkers),
                maxshape=(None, nwalkers),
                dtype=self.dtype,
            )

    def has_blobs(self):
        if self._has_blobs is None:
            raise RuntimeError('HDF5 backend not ready for use.')

        return self._has_blobs

    def get_value(self, name, flat=False, thin=1, discard=0):
        if not self.initialized:
            raise AttributeError(
                "You must run the sampler with "
                "'store == True' before accessing the "
                "results"
            )
        with self.open() as f:
            g = f[self.name]
            iteration = g.attrs["iteration"]
            if iteration <= 0:
                raise AttributeError(
                    "You must run the sampler with "
                    "'store == True' before accessing the "
                    "results"
                )

            assert g.attrs["has_blobs"] == self._has_blobs
            if name == "blobs" and not self._has_blobs:
                return None

            v = g[name][discard + thin - 1 : self.iteration : thin]
            if flat:
                s = list(v.shape[1:])
                s[0] = np.prod(v.shape[:2])
                return v.reshape(s)
            return v

    @property
    def shape(self):
        with self.open() as f:
            g = f[self.name]
            return g.attrs["nwalkers"], g.attrs["ndim"]

    @property
    def iteration(self):
        with self.open() as f:
            return f[self.name].attrs["iteration"]

    @property
    def accepted(self):
        with self.open() as f:
            return f[self.name]["accepted"][...]

    @property
    def random_state(self):
        with self.open() as f:
            elements = [
                v
                for k, v in sorted(f[self.name].attrs.items())
                if k.startswith("random_state_")
            ]
        #Inherited from EMCEE
        #pylint: disable=len-as-condition
        return elements if len(elements) else None
        #pylint: enable=len-as-condition
    #pylint: enable=missing-function-docstring

    def grow(self, ngrow, blobs):
        """Expand the storage space by some number of samples
        Args:
            ngrow (int): The number of steps to grow the chain.
            blobs: The current list of blobs. This is used to compute the
                dtype for the blobs array.
        """
        self._check_blobs(blobs)

        with self.open('r+') as f:
            g = f[self.name]
            ntot = g.attrs["iteration"] + ngrow
            g["chain"].resize(ntot, axis=0)
            g["log_prob"].resize(ntot, axis=0)
            if blobs is not None:
                assert self._has_blobs == g.attrs["has_blobs"]
                if not self._has_blobs:
                    nwalkers = g.attrs["nwalkers"]
                    dt = np.dtype((blobs[0].dtype, blobs[0].shape))
                    g.create_dataset(
                        "blobs",
                        (ntot, nwalkers),
                        maxshape=(None, nwalkers),
                        dtype=dt,
                    )
                else:
                    g["blobs"].resize(ntot, axis=0)
                self._has_blobs = True
                g.attrs["has_blobs"] = True
    #pylint: enable=invalid-name

    def save_step(self, state, accepted):
        """Save a step to the backend
        Args:
            state (State): The :class:`State` of the ensemble.
            accepted (ndarray): An array of boolean flags indicating whether
                or not the proposal for each walker was accepted.
        """

        retry_check = 10
        while retry_check > 0:
            try:
                self._check(state, accepted)
                retry_check = 0
            except BlockingIOError:
                retry_check -= 1
                if retry_check > 0:
                    logging.getLogger(__name__).error(
                        'Failed to _check step for saving. Retry in 1 min'
                    )
                    sleep(60)
                else:
                    raise

        with open(self.unsaved_steps_fname, 'ab') as unsaved_steps_file:
            pickle.dump((state, accepted), unsaved_steps_file)

        try:
            self._flush_unsaved_steps()
        except BlockingIOError:
            logging.getLogger(__name__).error(
                'Failed to save step to HDF5 file, will try again later'
            )
