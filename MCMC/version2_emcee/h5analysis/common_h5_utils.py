import numpy
import emcee
import sys
import scipy

from pathlib import Path
from directories import directories

home_dir=str(Path.home())
path=directories(home_dir)
sys.path.append(path.poet_path+'/PythonPackage')
sys.path.append(path.poet_path+'/scripts')

from orbital_evolution.transformations import phase_lag, lgQ

_blob_names = ['primary_mass',
               'secondary_mass', 
               'feh',
               'age',
               'eccentricity',
               'phase_lag_max',
               'alpha',
               'break_period']

_priors = {
    "primary_mass": {
        "min": 0.4,
        "max": 1.2
    },
    "secondary_mass": {
        "min": 0.4,
        "max": 1.2
    },
    "feh": {
        "min": 0.4,
        "max": 1.2
    },
    "age": {
        "min": 1e-2,
        "max": 10
    },
    "eccentricity": {
        "min": 0.0,
        "max": 0.8
    },
    "reference_lag": {
        "min": 2.984155182973038e-13,
        "max": 2.984155182973038e-06
    },
    "reference_lgQ": {
        "min": 5,
        "max": 12
    },
    "alpha": {
        "min": -5,
        "max": 5
    },
    "omega_break": {
        "min": numpy.log(2*numpy.pi/50),
        "max": numpy.log(2*numpy.pi/0.5)
    }
}


_tidal_periods = 10**numpy.linspace(numpy.log10(0.5), numpy.log10(50), 50)

_quantiles = [scipy.stats.norm.cdf(c) for c in [-2,-1,0,1,2]]

# def get_max_burn_in(kic):

#     burnins = []

#     with open(f'period_dependence/{kic}.txt','r') as f:
#         next(f)
#         for idx, lines in enumerate(f):
#             x = lines.split()
#             for i in [3,6,9,12]:
#                 burnins.append(int(x[i].split('/')[0]))
#     return max(burnins)

def get_max_burn_in(kic):

    burnins = numpy.empty((50,4), dtype=int)
    r_grid = numpy.empty((50,4), dtype=float)

    with open(f'period_dependence/{kic}.txt','r') as f:
        next(f)
        for idx, lines in enumerate(f):
            x = lines.split()
            burnins[idx] = numpy.array([int(x[i].split('/')[0]) for i in [3, 6, 9, 12]])
            r_grid[idx] = numpy.array([float(x[k]) for k in [2, 5, 8, 11]])

    max_burnin = max(burnins.flatten()[numpy.isfinite(r_grid.flatten())])
    return max_burnin

def best_constraint_period(kic):

    min_diff = 1e5
    with open(f'period_dependence/{kic}.txt','r') as f:
        next(f)
        for idx, lines in enumerate(f):
            x = lines.split()
            diff = float(x[10]) - float(x[1])
            period = float(x[0])
            if idx == 0 or diff < min_diff:
                min_period = period
                min_diff = diff
    
    return 10**min_period


def finite_posterior_samples(kic):

    old_h5filename = f'eightDh5files/system_{kic}.h5'
    new_h5filename = f'system_{kic}.h5'

    posterior_samples = dict()

    reader = emcee.backends.HDFBackend(new_h5filename, read_only=True)
    log_prob = reader.get_log_prob()
    mask = numpy.all(numpy.isfinite(log_prob), axis=1)
    blobs = reader.get_blobs()
    for name in _blob_names:
        posterior_samples[name] = blobs[name][mask]
    posterior_samples['log_prob'] = log_prob[mask]

    if numpy.all(mask):

        reader = emcee.backends.HDFBackend(old_h5filename, read_only=True)
        log_prob = reader.get_log_prob()
        mask = numpy.all(numpy.isfinite(log_prob), axis=1)
        blobs = reader.get_blobs()

        if numpy.any(mask):
            
            if numpy.all(mask): mask_idx = 0
            else:
                t = [i for i, x in enumerate(mask) if not x]
                mask_idx = max(t) + 1
            
            for name in _blob_names:

                if name == 'primary_mass':
                    blobs_val = blobs['m_sum']
                elif name == 'secondary_mass':
                    blobs_val = blobs['mass_ratio']
                elif name == 'feh':
                    blobs_val = blobs['metallicity']
                else:
                    blobs_val = blobs[name]
                
                posterior_samples[name] = numpy.concatenate((posterior_samples[name], blobs_val[mask_idx:]), axis=0)
            posterior_samples['log_prob'] = numpy.concatenate((posterior_samples['log_prob'], log_prob[mask_idx:]), axis=0)

    return posterior_samples

def converged_samples(kic, sample_dict):
    
    max_burn_in = get_max_burn_in(kic)
    for name in _blob_names:
        sample_dict[name] = sample_dict[name][max_burn_in:]
    sample_dict['log_prob'] = sample_dict['log_prob'][max_burn_in:]

def _power_law_funcion(delta0, 
                       alpha, 
                       omega_br, 
                       omega = 2*numpy.pi/20,
                       omega_min = 2*numpy.pi/50,
                       modify = False):

    if alpha > 0:
        if omega <= omega_min:
            delta = (omega_min/omega_br)**alpha
        elif omega > omega_min and omega <= omega_br:
            delta = (omega/omega_br)**alpha
        else:
            delta = 1
        delta *= delta0/((omega_min/omega_br)**alpha)
    else:
        if omega <= omega_br:
            delta = 1
        else:
            delta = (omega/omega_br)**alpha
        delta *= delta0

    return delta

def get_dissipation_samples(posterior_samples,
                            omega = 2*numpy.pi/20,
                            omega_min = 2*numpy.pi/50):

    sample_reference_lag = posterior_samples['reference_lag']
    samples_alpha = posterior_samples['alpha']
    samples_omega_br = posterior_samples['break_period']

    dissipation_samples = numpy.empty(numpy.shape(sample_reference_lag))
    lgQ_samples = numpy.empty(numpy.shape(sample_reference_lag))

    for i, row in enumerate(dissipation_samples):
        for j in range(len(row)):
            dissipation_samples[i,j] = _power_law_funcion(sample_reference_lag[i,j],
                                                          samples_alpha[i,j],
                                                          samples_omega_br[i,j],
                                                          omega = omega,
                                                          omega_min = omega_min)
            lgQ_samples[i,j] = lgQ(dissipation_samples[i,j])
    
    return (dissipation_samples, lgQ_samples)


