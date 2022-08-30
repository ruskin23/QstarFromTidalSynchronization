from itertools import chain
import emcee
from configargparse import ArgumentParser
import corner
import matplotlib.pyplot as plt
import numpy

#11704044

_quantities=['primary_mass',
             'secondary_mass',
             'feh',
             'age',
             'eccentricity',
             'phase_lag_max',
             'alpha',
             'omegeref']

def cmdline_args():

    p = ArgumentParser()
    p.add_argument('-s',
                    dest = 'chains_file')

    p.add_argument('-p',
                   dest = 'parameter')
    
    p.add_argument('--all',
                    action='store_true')

    p.add_argument('--steps',
                   action='store_true')
    
    return p.parse_args()



def get_samples(reader,
               parameter=None):
    
    samples=reader.get_blobs()


def corner_plot(reader,
                samples=None):

    samples = reader.get_chain(flat=True)
    log_prob_samples = reader.get_log_prob(flat=True)
    log_prior_samples = reader.get_blobs(flat=True) 

    print(numpy.shape(samples),
          numpy.shape(log_prob_samples),
          numpy.shape(log_prior_samples))

    # all_samples = numpy.concatenate((samples, log_prob_samples[:, None], log_prior_samples[:, None]), axis=1)

    corner.corner(samples)
    plt.show()




if __name__=='__main__':

    args = cmdline_args()
    
    filename = args.chains_file
    parameter = args.parameter

    if parameter is not None and parameter not in _quantities:
        print('parameter defined not in predefined quantities: {}'.format(_quantities))

    reader = emcee.backends.HDFBackend(filename, read_only=True)

    if args.steps:
        system_kic = filename.split('.')[0].split('_')[1]
        try:
            print('Total Steps taken for system {} = {}'.format(system_kic, len(reader.get_log_prob(flat=True))//64))
        except:
            print('error raised for system {}'.format(system_kic))


    corner_plot(reader)

