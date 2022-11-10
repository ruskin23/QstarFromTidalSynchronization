import numpy
import common_h5_utils
import corner
import matplotlib.pyplot as plt
import json
import sys

from pathlib import Path
from directories import directories

home_dir=str(Path.home())
path=directories(home_dir)
sys.path.append(path.poet_path+'/PythonPackage')
sys.path.append(path.poet_path+'/scripts')

from orbital_evolution.transformations import phase_lag, lgQ
from discrete_utils import discrete_stats



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
                            
    plt.savefig(f'plots/corner_{kic}.png')
    plt.close()


if __name__ == '__main__':    

    with open('convergence.json', 'r') as f:
        convergence_dict = json.load(f)

    for kic in convergence_dict.keys():
        if convergence_dict[kic]['converged'] == 'True':

            posterior_samples = common_h5_utils.finite_posterior_samples(kic)
            common_h5_utils.converged_samples(kic, posterior_samples)

            posterior_samples['reference_lag'] = posterior_samples['phase_lag_max']
            del posterior_samples['phase_lag_max']

            pdf = discrete_stats(posterior_samples)

            delta_samples, lgQ_samples = common_h5_utils.get_dissipation_samples(posterior_samples)
            

    # data = []
    # for keys, values in posterior_samples.items():
    #     print(len(values))
    #     print(numpy.isnan(numpy.sum(values)))
    # data = numpy.vstack(data).T
    # figure = corner.corner(data)
    # plt.savefig('temp2.png')