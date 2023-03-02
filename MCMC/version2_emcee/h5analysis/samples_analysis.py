import numpy
import sys
import pickle
import json
import os
import logging

from multiprocessing import Pool

import argparse

import seaborn as sns
import matplotlib.pyplot as plt
import corner

from pathlib import Path
from result_analysis_directories import directories

home_dir=str(Path.home())
path=directories(home_dir)
sys.path.append(path.poet_path+'/PythonPackage')
sys.path.append(path.poet_path+'/scripts')
sys.path.append(path.mcmc_directory)
sys.path.append('/home/ruskin/projects/')

_working_directory = path.mcmc_directory + '/h5analysis'

from orbital_evolution.transformations import phase_lag, lgQ
import utils
import common_h5_utils
from sampler_logger import setup_logging
from general_purpose_python_modules.emcee_quantile_convergence import *



_joint_params = ['alpha', 'omega_break', 'reference_lgQ']

_quantitites = ['primary_mass',
               'secondary_mass',
               'feh',
               'age',
               'eccentricity',
               'reference_lag',
               'reference_lgQ',
               'alpha',
               'break_period']
tex_fonts = {
    # Use LaTeX to write all text
    "text.usetex": True,
    "font.family": "serif",
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 15,
    # "font.size": 15,
    "font.weight" : "bold",
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 20,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12
}
plt.rcParams.update(tex_fonts)

_logger=logging.getLogger(__name__)


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

    parser.add_argument('--npriors',
                    type=int,
                    default=1000,
                    help='dimensinon of prior grid')

    parser.add_argument('--logging_path',
                    default=path.scratch_directory + '/sampling_output/sampler',
                    help='output directory')

    parser.add_argument('--std_out_err_path',
                    default=path.scratch_directory + '/sampling_output/sampler',
                    help='output directory')

    parser.add_argument('--sampler',
                        action='store_true',
                        default = False,
                        help='get samples or nah')


    parser.add_argument('--pickle',
                         action='store_true',
                         default = False,
                         help='use this if there is a pickle for posterior available')

    return parser.parse_args()

def corner_plot(kic, posterior_samples):

    data = numpy.transpose(
                            numpy.vstack(
                                            (
                                                (
                                                    posterior_samples[kic]['age']['samples'].flatten(),
                                                    posterior_samples[kic]['reference_lgQ']['samples'].flatten(),
                                                    posterior_samples[kic]['alpha']['samples'].flatten(),
                                                    2*numpy.pi/posterior_samples[kic]['omega_break']['samples'].flatten()
                                                )
                                            )
                                        )
                            )

    figure = corner.corner(data,
                            labels=[
                            r"age",
                            r"$\log_{10}{Q_{\ast,0}^{\prime}}$",
                            r"$\alpha$",
                            r"$\log_{10}{P_{br}}$"
                            ],
                            quantiles=[0.16, 0.5, 0.84],
                            show_titles=True,
                            title_kwargs={'fontsize': 12},
                            )

    # plt.savefig(_working_directory + f'/plots/corner_{kic}.png')
    plt.savefig('corner8543278.pdf')
    plt.close()

    if 'logQ_best' in posterior_samples.keys():
        sns.jointplot(x=posterior_samples["alpha"], y=posterior_samples["logQ_best"], kind='hex')
        plt.savefig(_working_directory + f'/plots/2D_alphaQ_{kic}.png')
        plt.close()

        sns.jointplot(x=posterior_samples["omega_break"], y=posterior_samples["logQ_best"], kind='hex')
        plt.savefig(_working_directory + f'/plots/2D_breakQ_{kic}.png')
        plt.close()


def create_posterior_dataset(parse_args):

    if parse_args.pickle:
        with open(_working_directory + '/posterior_dataset.pickle', 'rb') as f:
            return pickle.load(f)

    with open(_working_directory + '/convergence.json', 'r') as f:
        convergence_dict = json.load(f)

    prior_limits = common_h5_utils._priors

    distribution_dict = {}

    _quantitites[_quantitites.index('break_period')] = 'omega_break'

    # converged = ['4380283', '11147276', '10091257', '6525196', '5022440', '5802470', '4815612', '3241344', '10031409', '7838639', '8957954', '4285087', '11233911', '12004679', '9971475', '7362852', '8543278', '6927629', '8364119', '6949550', '3348093', '9892471', '4839180', '5652260', '11232745', '8984706', '5393558', '9353182', '7336754', '6029130', '11704044', '4678171', '7987749', '10330495', '10992733', '10711913', '10215422', '4773155']

    for kic in convergence_dict.keys():
        if convergence_dict[kic]['converged'] == 'True':
    # for kic in converged:
    #     if True:
            distribution_dict[kic] = dict()

            posterior_samples = common_h5_utils.finite_posterior_samples(kic)
            common_h5_utils.converged_samples(kic, posterior_samples)
            print(kic)
            print(len(posterior_samples['phase_lag_max']))

            posterior_samples['reference_lag'] = posterior_samples['phase_lag_max']
            posterior_samples['omega_break'] = posterior_samples['break_period']
            posterior_samples['reference_lgQ'] = numpy.array([lgQ(lag) for lag in posterior_samples['reference_lag'].flatten()]).reshape((numpy.shape(posterior_samples['reference_lag'])))
            del posterior_samples['phase_lag_max']
            del posterior_samples['break_period']

            distribution_dict[kic]['n_samples'] = len(posterior_samples['reference_lag'].flatten())

            for names in _quantitites:
                distribution_dict[kic][names] = dict()
                distribution_dict[kic][names]['samples'] = posterior_samples[names].flatten()
                distribution_dict[kic][names]['prior'] = numpy.linspace(prior_limits[names]['min'],
                                                                         prior_limits[names]['max'],
                                                                         parse_args.npriors)
                if names == 'omega_break':
                    distribution_dict[kic][names]['prior'] = numpy.exp(distribution_dict[kic][names]['prior'])
                distribution_dict[kic][names]['bandwidth'] = utils._get_kernel_bandwidth(posterior_samples[names].flatten())


    with open(_working_directory + '/posterior_dataset.pickle', 'wb') as f:
        pickle.dump(distribution_dict, f)
    return distribution_dict


class joint_distribution():

    def __init__(self,
                 posterior_dataset,
                 n_prior = 1000):

        self._logger = logging.getLogger(__name__)
        self.posterior_dataset = posterior_dataset.copy()
        self.n_prior = n_prior

    def update_weights(self, weights, sampled):

        for kic in self.posterior_dataset.keys():
            weights[kic] *= self.posterior_dataset[kic]['updated_weight_func'].get_norm(sampled)


    def get_sample(self, unit_vector):

        _logger.info(f'\n\nSAMPLING RANDOM VALUE AT UNIT VECTOR: {unit_vector}')
        unit_vector_iter = iter(unit_vector)
        sampled_tup = []

        weights = {}
        for kic in self.posterior_dataset.keys():
            weights[kic] = numpy.ones(self.posterior_dataset[kic]['n_samples'])

        for param in _joint_params:

            self.param_pdf = numpy.ones(self.n_prior, dtype=float)

            _logger.info(f'\nCalculating for {param}')
            for kic, sample_set in self.posterior_dataset.items():

                _logger.info(f'\nCalculating for {kic}')
                #1 create a joint pdf at some prior value array for each parameter
                distribution_func = utils.DiscreteSampling(sample_set[param]['samples'], sample_set[param]['bandwidth'])

                _logger.info('Assignging weights')
                distribution_func.set_weights(weights[kic])
                self.posterior_dataset[kic]['updated_weight_func'] = distribution_func

                _logger.info('Updating pdf grid')
                for i, prior_val in enumerate(sample_set[param]['prior']):
                    self.param_pdf[i] *= distribution_func._pdf(prior_val)

            #2 normalize pdf
            _logger.info('Normalizing join pdf')
            normalization = utils._normalization_constant(sample_set[param]['prior'], self.param_pdf, min(sample_set[param]['prior']), max(sample_set[param]['prior']))
            self.param_pdf /= normalization

            #3 sample a value corresponding to unit cube from the discrete joint pdf
            _logger.info('Sampling at random unit value')
            discrete_distribution = utils.DiscreteDistributions(self.param_pdf, sample_set[param]['prior'])
            param_sampled = discrete_distribution._ppf(next(unit_vector_iter))
            sampled_tup.append(param_sampled)

            #4Update weigths for individual kic given the sampled parameter for the next param
            _logger.info('Updating Weights')
            self.update_weights(weights, param_sampled)

        _logger.info('Sampling Complete')
        return sampled_tup

def sample_params(parse_args):

    posterior_dataset = create_posterior_dataset(parse_args)

    unit_vector = [numpy.random.rand(len(_joint_params)) for i in range(parse_args.nsamples)]

    joint = joint_distribution(posterior_dataset, n_prior=parse_args.npriors)

    with Pool(processes=parse_args.nprocs,
              initializer=setup_logging,
              initargs=[parse_args],
              maxtasksperchild=1
              ) as pool:

              sampled = pool.map(joint.get_sample, unit_vector)

    with open(_working_directory + f'/sampled_params{parse_args.nsamples}.pickle', 'wb') as f:
        pickle.dump(sampled, f)

def plot_parameter_corner(posterior_samples, pasrse_args):

    with open(_working_directory + '/convergence.json', 'r') as f:
        convergence_dict = json.load(f)

    for kic in convergence_dict.keys():
        if convergence_dict[kic]['converged'] == 'True' and kic != '5393558':

            print(f'Plotting for kic {kic}')
            alpha = posterior_samples[kic]['alpha']['samples'].flatten()
            omega_break = 2*numpy.pi/posterior_samples[kic]['omega_break']['samples'].flatten()
            omega_break = numpy.log10(omega_break)
            reference_lag = posterior_samples[kic]['reference_lag']['samples'].flatten()
            lgQ_reference = numpy.array([lgQ(lag) for lag in reference_lag])
            data = numpy.transpose(
                                    numpy.vstack(
                                                    (
                                                        (
                                                            lgQ_reference,
                                                            alpha,
                                                            omega_break
                                                        )
                                                    )
                                                )
                                    )
            figure = corner.corner(
                                    data,
                                    axes_scale='log, linear, log',
                                    labels=[
                                        r"$\log_{10}{Q_{\ast,0}^{\prime}}$",
                                        r"$\alpha$",
                                        r"$\log_{10}{P_{br}}$"
                                        ],
                                        quantiles=[0.16, 0.5, 0.84],
                                        show_titles=True,
                                        title_kwargs={'fontsize': 15},
                                        )

            plt.savefig(f'plots/alpha_break_{kic}.png')
            plt.savefig(f'plots/alpha_break_{kic}.pdf')
            plt.close()


    with open(f'sampled_params{parse_args.nsamples}.pickle', 'rb') as f:
        sampled = pickle.load(f)

    data = []
    for val in sampled:
        data.append((val[2], val[0], 2*numpy.pi/val[1]))

    figure = corner.corner(
                            data,
                            axes_scale='log, linear, log',
                            labels=[
                                r"$\log_{10}{Q_{\ast,0}^{\prime}}$",
                                r"$\alpha$",
                                r"$\log_{10}{P_{br}}$"
                                ],
                                quantiles=[0.16, 0.5, 0.84],
                                show_titles=True,
                                title_kwargs={'fontsize': 15, 'weight': 'bold'},
                                )

    plt.savefig(f'plots/alpha_break_combined.png')
    plt.savefig(f'plots/alpha_break_combined.pdf')
    plt.tight_layout()
    plt.close()


def rl_convergence_test(kic):

    _logger = logging.getLogger(__name__)

    _logger.info(f'Testing for KIC {kic}')
    posterior_samples = common_h5_utils.finite_posterior_samples(kic)

    posterior_samples['reference_lag'] = posterior_samples['phase_lag_max']
    del posterior_samples['phase_lag_max']
    
    stats_grid = numpy.ones((len(common_h5_utils._quantiles), len(common_h5_utils._tidal_periods), 4), dtype=float)

    max_step = 0
    max_burning = 0

    for q_idx, q in enumerate(common_h5_utils._quantiles):
        _logger.info(f'\nCalculating for quantile {q}')

        temp_grid = numpy.zeros((len(common_h5_utils._tidal_periods), 4,))

        for p_idx, period in enumerate(common_h5_utils._tidal_periods):

            _logger.info(f'At Period = {period} days')
            omega = 2*numpy.pi/period
            _, lgQ_samples = common_h5_utils.get_dissipation_samples(posterior_samples, omega=omega)

            if max_step == 0: max_step = len(lgQ_samples)

            quantile_val, _, std, thinning, burnin =  find_emcee_quantiles(lgQ_samples, q , 0.001, 1000)

            max_burning = max(max_burning, burnin)
            
            temp_grid[p_idx] = numpy.array([quantile_val, std, thinning, burnin])
        
        stats_grid[q_idx] = temp_grid

    if max_burning > max_step: _logger.info(f'{kic} NOT_CONVERGED')
    elif max_step - max_burning < 100: _logger.info(f'{kic} NOT_CONVERGED LOW_STEPS')
    else: _logger.info(f'{kic} CONVERGED')

    with open(f'rl_stats_dump/{kic}.pickle', 'wb') as f:
        pickle.dump(stats_grid, f)





def test_convergence(parse_args):

    with open('convergence.json', 'r') as f:
        system_test = json.load(f)
    
    # kic_list = []
    # for kic, stat in system_test.items():
    #     if stat['converged'] == 'False':kic_list.append(kic)

    kic_list = ['5022440', '4352168', '11616200', '9656543']

    with Pool(processes=parse_args.nprocs,
              initializer=setup_logging,
              initargs=[parse_args],
              maxtasksperchild=1
              ) as pool:

              pool.map(rl_convergence_test, kic_list)

def convergence_stats(kic):

    with open(f'rl_stats_dump/{kic}.pickle', 'rb') as f:
        rl_stats = pickle.load(f)

    max_steps = len(common_h5_utils.finite_posterior_samples(kic)['log_prob'])

    pdf_grid = numpy.zeros((len(common_h5_utils._tidal_periods), 1000), dtype=float)
    cdf_grid = numpy.zeros((len(common_h5_utils._tidal_periods), 1000), dtype=float)

    
    

    




if __name__ == '__main__':

    parse_args = cmd_parser()
    # posterior = create_posterior_dataset(parse_args)
    # corner_plot('8543278', posterior)
    # print(posterior['5393558']['reference_lgQ'])
    if parse_args.sampler: sample_params(parse_args)
    # plot_parameter_corner(posterior, parse_args)

    # test_convergence = 

    # rl_convergence_test('10711913')
    test_convergence(parse_args)




