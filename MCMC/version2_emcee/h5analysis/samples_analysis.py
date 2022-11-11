import numpy
import sys
import pickle
import json

from multiprocessing import Pool

import argparse
import logging

import seaborn as sns
import matplotlib.pyplot as plt
import corner

from pathlib import Path
from directories import directories

home_dir=str(Path.home())
path=directories(home_dir)
sys.path.append(path.poet_path+'/PythonPackage')
sys.path.append(path.poet_path+'/scripts')
sys.path.append(path.mcmc_directory)
_working_directory = path.mcmc_directory + '/h5analysis'

from orbital_evolution.transformations import phase_lag, lgQ
import utils
from Logger import setup_process 
import common_h5_utils


_joint_params = ['alpha', 'break_period']

_blob_names = ['primary_mass',
               'secondary_mass', 
               'feh',
               'age',
               'eccentricity',
               'reference_lag',
               'alpha',
               'break_period']


def cmd_parser():

    parser = argparse.ArgumentParser()
    
    parser.add_argument('--nprocs', 
                        default=16,
                        type=int,
                        help='number of processes')

    parser.add_argument('--nsamples',
                        type=int,
                        default=256,
                        help='number of samples')

    parser.add_argument('--fname_datetime_format', default='%Y%m%d%H%M%S')
    parser.add_argument('--logging_datetime_format', default=None)
    parser.add_argument('--std_out_err_fname', default='%(system)s_%(now)s_%(pid)d.err')
    parser.add_argument('--logging_message_format', default=None)
    parser.add_argument('--logging_fname', default='%(system)s_%(now)s_%(pid)d.log')

    return parser.parse_args()

def corner_plot(kic, posterior_samples):

    data = numpy.transpose(
                            numpy.vstack(
                                            (
                                                (
                                                    numpy.array([lgQ(lag) for lag in posterior_samples['reference_lag'].flatten()]),
                                                    posterior_samples['alpha'].flatten(), 
                                                    numpy.log10(2*numpy.pi/posterior_samples['break_period'].flatten())
                                                )
                                            )
                                        )
                            )

    figure = corner.corner(data,
                            labels=[
                            r"$\log_{10}{Q_{\ast}^{\prime}}$",
                            r"$\alpha$",
                            r"$\log_{10}{P_{br}}$"
                            ],
                            quantiles=[0.16, 0.5, 0.84],
                            show_titles=True,
                            title_kwargs={'fontsize': 12},
                            )

    plt.savefig(_working_directory + f'/plots/corner_{kic}.png')
    plt.close()

    if 'logQ_best' in posterior_samples.keys():
        sns.jointplot(x=posterior_samples["alpha"], y=posterior_samples["logQ_best"], kind='hex')
        plt.savefig(_working_directory + f'/plots/2D_alphaQ_{kic}.png')
        plt.close()

        sns.jointplot(x=posterior_samples["break_period"], y=posterior_samples["logQ_best"], kind='hex')
        plt.savefig(_working_directory + f'/plots/2D_breakQ_{kic}.png')
        plt.close()


class joint_distribution():

    def __init__(self, 
                 posterior_dataset,
                 n_prior = 1000):

        self.posterior_dataset = posterior_dataset
        self.n_prior = n_prior
        self.param_pdf = numpy.ones(n_prior, dtype=float)
        
    def update_weights(self, weights, sampled):

        for kic in self.posterior_dataset.keys():
            weights[kic] *= self.posterior_dataset[kic]['updated_weight_func'].get_norm(sampled)


    def get_sample(self, unit_vector):

        _logger = logging.getLogger(__name__)

        unit_vector_iter = iter(unit_vector)
        sampled_set = []

        weights = {}
        for kic in self.posterior_dataset.keys():
            weights[kic] = numpy.ones(self.posterior_dataset[kic]['n_samples'])

        for param in _joint_params:
            _logger.info(f'\nCalculating for param = {param}')

            for kic, sample_set in self.posterior_dataset.items():
                _logger.info(f'\nfor kic = {kic}')
                #1 create a joint pdf at some prior value array for each parameter
                distribution_func = utils.DiscreteSampling(sample_set[param]['samples'], sample_set[param]['bandwidth'])
                distribution_func.set_weights(weights[kic])
                self.posterior_dataset[kic]['updated_weight_func'] = distribution_func
                self.param_pdf *= distribution_func._pdf(sample_set[param]['prior'])

            #2 normalize pdf
            normalization = utils._normalization_constant(sample_set[param]['prior'], self.param_pdf, min(sample_set[param]['prior']), max(sample_set[param]['prior']))
            self.param_pdf /= normalization

            #3 sample a value corresponding to unit cube from the discrete joint pdf
            discrete_distribution = utils.DiscreteDistributions(self.param_pdf, sample_set[param]['prior'])
            param_sampled = discrete_distribution._ppf(next(unit_vector_iter))
            sampled_set.append(param_sampled)

            #4Update weigths for individual kic given the sampled parameter for the next param
            self.update_weights(weights, param_sampled)

        return sampled_set

def create_posterior_dataset(nsamples):

    _logger = logging.getLogger(__name__)
    with open(_working_directory + '/convergence.json', 'r') as f:
        convergence_dict = json.load(f)

    with open(_working_directory + '/priors.json', 'r') as f:
        prior_limits = json.load(f)

    distribution_dict = {}

    _logger.info('Creating Dataset')
    for kic in convergence_dict.keys():
        if convergence_dict[kic]['converged'] == 'True':
            
            _logger.info(f'\n\n{kic}')

            distribution_dict[kic] = dict()

            posterior_samples = common_h5_utils.finite_posterior_samples(kic)
            common_h5_utils.converged_samples(kic, posterior_samples)

            posterior_samples['reference_lag'] = posterior_samples['phase_lag_max']
            del posterior_samples['phase_lag_max']
            distribution_dict[kic]['n_samples'] = len(posterior_samples['reference_lag'].flatten())

            for names in _blob_names:
                distribution_dict[kic][names] = dict()
                distribution_dict[kic][names]['samples'] = posterior_samples[names].flatten()
                distribution_dict[kic][names]['prior'] = numpy.linspace(prior_limits[names]['min'], 
                                                                         prior_limits[names]['max'],
                                                                         nsamples)
                distribution_dict[kic][names]['bandwidth'] = utils._get_kernel_bandwidth(posterior_samples[names].flatten())

    return distribution_dict

def sample_params(parse_args):

    setup_process(parse_args)

    posterior_dataset = create_posterior_dataset(parse_args.nsamples)

    unit_vector = [numpy.random.rand(2) for i in range(parse_args.nsamples)]

    joint = joint_distribution(posterior_dataset, n_prior=parse_args.nsamples)

    with Pool(processes=parse_args.nprocs,
              initializer=setup_process,
              initargs=[parse_args],
              maxtasksperchild=1
              ) as pool:
              
              sampled = pool.map(joint.get_sample, unit_vector)

    with open(_working_directory + '/sampled_params.pickle', 'wb') as f:
        pickle.dump(sampled, f)


if __name__ == '__main__':    

    parse_args = cmd_parser()
    sample_params(parse_args)