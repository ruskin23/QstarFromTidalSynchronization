from KDEpy import FFTKDE
import emcee
import sys
from configargparse import ArgumentParser
import corner
import matplotlib.pyplot as plt
import numpy
import h5py
import glob

from pathlib import Path
from directories import directories


home_dir=str(Path.home())
path=directories(home_dir)
sys.path.append(path.poet_path+'/PythonPackage')
sys.path.append(path.poet_path+'/scripts')
sys.path.append('/home/ruskin/projects/')
sys.path.append('/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/version2_emcee/')

from orbital_evolution.transformations import phase_lag, lgQ
from general_purpose_python_modules.emcee_quantile_convergence import *
from utils import *
import scipy
import scipy.stats


def save_systems():

    rerun_systems = []
    remaining_systems = []
    skipped_steps_systems = []

    saved_systems_files = ['rerun_systems.out', 'remaining_systems.out', 'skipped_steps.out']
    for sf in saved_systems_files:
        with open(sf,'w') as f:
            pass


    filenames = glob.glob('*.h5')
    for files in filenames:
        system_kic = files.split('.')[0].split('_')[1]

        reader = emcee.backends.HDFBackend(files, read_only=True)
        log_prob = reader.get_log_prob()
        
        data = h5py.File(files, 'r')
        all_log_probs = data['/mcmc/log_prob']
        check = numpy.isfinite(all_log_probs).all()

        if not check:
            total_steps = len(log_prob)
            for step in reversed(range(total_steps)):
                if -numpy.inf in log_prob[step]:
                    break            
            # extra_check = numpy.isfinite(all_log_probs[step:total_steps].all())
            print('For system {} step = {} total_steps = {}'.format(system_kic, step+1, total_steps))
            if step == total_steps-1:
                with open('rerun_systems.out','a') as f:
                    f.write(system_kic+'\n')
                rerun_systems.append(system_kic)
            else:
                with open('skipped_steps.out','a') as f:
                    f.write(system_kic+'\t'+repr(step)+'\t'+repr(total_steps)+'\n')
                skipped_steps_systems.append(system_kic)

    catalog = '/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/version2_emcee/catalog/filtering/nominal_value_catalog_Iconv_cutoff.txt'

    with open(catalog,'r') as f:
        next(f)
        for lines in f:
            x = lines.split()
            system_kic = x[1]
            if system_kic not in rerun_systems and system_kic not in skipped_steps_systems:
                with open('remaining_systems.out','a') as f:
                    f.write(system_kic+'\n')
                remaining_systems.append(system_kic)

    print('\nRerun Systems:')
    print(' '.join(rerun_systems))

    print('\nSkipped Steps Systems')
    print(' '.join(skipped_steps_systems))

    print('\nRemaining Systems:')
    print(' '.join(remaining_systems))

    print('\n',len(rerun_systems),len(skipped_steps_systems),len(remaining_systems))


def get_possible_last_step_from_oldfiles():

    with open('continued_systems_recheck.txt','w') as f:
        f.write('kic'+'\t'+
                'inf_present'+'\t'+
                'finite_step'+'\t'+
                'total_steps'+'\n')
    continued_systems = []
    with open('remaining_systems.out','r') as f:
        for lines in f:
            x=lines.split()
            continued_systems.append((x[0]))

    filenames = glob.glob('eightDh5files/*.h5')

    for f_name in filenames:
        system_kic = f_name.split('.')[0].split('_')[1]
        if system_kic in continued_systems:

            data = h5py.File(f_name, 'r')
            all_log_probs = data['/mcmc/log_prob']
            check = numpy.isfinite(all_log_probs).all()

            reader = emcee.backends.HDFBackend(f_name, read_only=True)
            log_prob = reader.get_log_prob()

            total_steps = len(log_prob)

            if not check:

                for step in reversed(range(total_steps)):
                    if -numpy.inf in log_prob[step]:
                        break

            else:
                step = 0
            
            with open('continued_systems_recheck.txt','a') as f:
                f.write(system_kic + '\t' +
                        repr(check) + '\t' +
                        repr(step) + '\t' +
                        repr(total_steps) + '\n')




def corner_plot(system_KIC,
                samples):

    corner.corner(samples,
                  labels=[
                        r"m$_p$",
                        r"m$_s$",
                        r"feh",
                        r"age",
                        r"$e$",
                        r"$\log_{10}(\Delta)$",
                        r"$\alpha$",
                        r"$\omega_{br}$",
                        ],
                  quantiles=[0.1, 0.5, 0.9],
                  show_titles=True,
                  title_kwargs={"fontsize": 10},
                )

    plt.show()
    plt.savefig('corner_plots/{}.png'.format(system_KIC))


def rearrange_prior_samples(samples,
                            truncate=True,
                            dimensions=9):

    samples = samples.flatten()
    rearranged_samples = [[samples[k][i] for i in range(dimensions)] for k in range(len(samples))]
    if truncate: rearranged_samples = numpy.delete(rearranged_samples, 5, axis=1)

    return rearranged_samples

def corrected_samples(system_filename):

    system_kic =  system_filename.split('.')[0].split('_')[1]
    reader = emcee.backends.HDFBackend(system_filename,read_only=True)

    include_old = False
    with open('fix_steps.txt','r') as f:

        for _ in range(3):
            next(f)

        for lines in f:
            x = lines.split()
            if x[0] == system_kic:
                if x[1] == 'True':
                    prior_samples = reader.get_blobs()
                elif x[2] == 'True':
                    skip_steps = int(x[3])
                    prior_samples = reader.get_blobs()[skip_steps:]
                else:
                    include_old = True
                    prior_samples = reader.get_blobs()
                    reader_old = emcee.backends.HDFBackend('eightDh5files/'+system_filename, read_only=True)
                    if x[4] == 'False':
                        skip_steps = int(x[5])
                        prior_samples_old = reader_old.get_blobs()[skip_steps:]
                    elif x[4] == 'True':
                        prior_samples_old = reader_old.get_blobs()
                    else:
                        raise ValueError('None of the conditions satisfied')

    if include_old:
        samples_old = rearrange_prior_samples(prior_samples_old, truncate=False, dimensions=8)
        sample_new = rearrange_prior_samples(prior_samples)
        samples = numpy.concatenate((samples_old, sample_new), axis=0)
    else:
        samples = rearrange_prior_samples(prior_samples)

    return samples    

def cmdline_args():

    p = ArgumentParser()
    p.add_argument('-s',
                    dest = 'chains_file',
                    default=None)

    p.add_argument('-p',
                   dest = 'parameter')

    p.add_argument('-c',
                    dest='convergence_filename',
                    default=None)

    p.add_argument('--all',
                    action='store_true',
                    default=False)
    
    p.add_argument('--expand',
                    action='store_true',
                    default=False)

    p.add_argument('--steps',
                   action='store_true',
                   default=False)

    return p.parse_args()

def get_dissipation_samples(samples, tidal_period=10):

    q_samples = []
    for s in samples:
        delta0 = s[5]
        alpha = s[6]
        omega_br = s[7]
        omega_min = 2*numpy.pi/50
        omega = 2*numpy.pi/tidal_period
        
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
        
        q_samples.append(lgQ(delta))

    if numpy.nan in q_samples:        
        print('NaN found in Q_SAMPLES!!!!!!!!!')
    q_samples = numpy.reshape(q_samples,(len(q_samples)//64,64))

    return q_samples

def emcee_convergence_test(lgQ_samples,
                           save_filename=None):

    total_steps = numpy.shape(lgQ_samples)[0]
    _quantities=[scipy.stats.norm.cdf(c) for c in [-2,-1,1,2]]
    rl_stat_list = []
    for qunatile in _quantities:
        print('Calculating for quantile {}'.format(qunatile))
        rl_stats = find_emcee_quantiles(lgQ_samples,qunatile,0.001,1000)
    
        qunatile_val = rl_stats[0]
        r = rl_stats[2]
        thin = rl_stats[3]
        burnin = rl_stats[4]
        if save_filename is not None:
            with open(save_filename,'a') as f:
                f.write('\t' + 
                        repr(qunatile) + '\t' +
                        repr(r) + '\t' +
                        repr(qunatile_val) + '\t' +
                        repr(burnin) + '\t' +
                        repr(total_steps) + '\t' +
                        repr(thin)+'\n')
        else: 
            print(rl_stats)
            rl_stat_list.append(rl_stats)

    return rl_stat_list

def tidal_period_dependence(system_kic, prior_samples):    
    
    quantiles = [scipy.stats.norm.cdf(c) for c in [-2,-1,1,2]]
    f = open('period_dependence/'+system_kic+'.txt', 'w', buffering=1)
    f.write('period'+'\t'+
            '-2sigma({})'.format(repr(quantiles[0]))+'\t'+
            'r'+'\t'+
            'burnin/total_steps'+'\t'+
            '-1sigma({})'.format(repr(quantiles[1]))+'\t'+
            'r'+'\t'+
            'burnin/total_steps'+'\t'+
            '1sigma({})'.format(repr(quantiles[2]))+'\t'+
            'r'+'\t'+
            'burnin/total_steps'+'\t'+
            '2sigma({})'.format(repr(quantiles[3]))+'\t'+
            'r'+'\t'+
            'burnin/total_steps'+'\n')

    period_limts = [numpy.log10(0.5),numpy.log10(50),50]
    tidal_period = numpy.linspace(float(period_limts[0]), float(period_limts[1]), int(period_limts[2]))
    data = numpy.empty((50,5), float)
    data_std = numpy.empty((50,4), float)
    for i, p in enumerate(tidal_period):
        lgQ_samples = get_dissipation_samples(prior_samples, tidal_period=10**p)
        rl_stats = emcee_convergence_test(lgQ_samples)

        line = '\t'.join(['\t'.join([repr(s[0]), repr(s[2]), '{}/{}'.format(repr(s[4]),repr(len(lgQ_samples)))]) for s in rl_stats])
        f.write(repr(p) + '\t' + line + '\n')

        data[i] = numpy.array([p] + [rl_stats[i][0] for i in range(4)])
        data_std[i] = numpy.array([rl_stats[i][2] for i in range(4)])
        
    f.close()

    #Plot converged systems
    if numpy.isnan(numpy.sum(data_std)):
        print('System Not Converged Yet')
        return
    else:
        with open('period_dependence/new_stop_systems.txt','a') as f:
            f.write('{}\n'.format(system_kic))
    data = numpy.transpose(data)
    for k in range(4):
        plt.plot(data[0,:i],data[k+1,:i])
    plt.savefig('period_dependence/plot_'+system_kic+'.png')


def plot_converged(system_kic=None):
    #Plot quantiles of Q vs Period all the converged systems or a specific system

    if system_kic is None:
        with open('period_dependence/new_stop_systems.txt','r') as f:
            system_kic = f.read().split('\n')
        if '' in system_kic:
            system_kic.remove('')

    for S in system_kic:
        data = numpy.empty((50,5), float)
        with open(f'period_dependence/{S}.txt','r') as f:
            next(f)
            for i,lines in enumerate(f):
                x = lines.split()
                data[i] = numpy.array([float(x[k]) for k in [0,1,4,7,10]])
    
        data = numpy.transpose(data)
        for k in range(4):
            plt.plot(data[0,:i],data[k+1,:i])
        plt.savefig('period_dependence/plot_'+S+'.png')
        plt.close()

def plot_max_quantile_difference():
    
    with open('period_dependence/new_stop_systems.txt','r') as f:
        converged_systems = f.read().split('\n')
    converged_systems = [sk for sk in converged_systems if sk != '']
    period_limts = numpy.linspace(numpy.log10(0.5),numpy.log10(50),50)
    
    q_diff_dict = dict()
    for system_kic in converged_systems:
        Q_diff = []
        with open(f'period_dependence/{system_kic}.txt', 'r') as f:
            next(f)
            for lines in f:
                x = lines.split()
                Q_diff.append(float(x[10]) - float(x[1]))
        q_diff_dict[system_kic] = Q_diff


def get_pdf(kic):

    burnins = []
    with open('period_dependence/{}.txt'.format(kic),'r') as f:
        next(f)
        for lines in f:
            x = lines.split()
            for i in [3,6,9,12]:
                burnins.append(int(x[i].split('/')[0]))
                max_step = int(x[i].split('/')[1])

    max_burnins = max(burnins)


    prior_samples = corrected_samples(f'system_{kic}.h5')
    periods = numpy.linspace(numpy.log10(0.5), numpy.log10(50), 50)
    lgQ_array = numpy.linspace(5,20,250)

    Q_diff = 0.2
    lgQ_sigma_grid = numpy.empty((50,5),float)
    with open('period_dependence/{}.txt'.format(kic),'r') as f:
        next(f)
        for i,lines in enumerate(f):
            x = lines.split()
            lgQ_sigma_grid[i] = numpy.array([[float(x[k]) for k in [0,1,4,7,10]]])
    lgQ_sigma_grid = numpy.transpose(lgQ_sigma_grid)
    lgQ_2sigma_upper = lgQ_sigma_grid[-1]
    lgQ_2sigma_lower = lgQ_sigma_grid[1]
    minQ = lgQ_2sigma_upper[numpy.argmin(lgQ_2sigma_upper - lgQ_2sigma_lower)]
    
    reduced_periods = lgQ_sigma_grid[0][lgQ_2sigma_upper < (minQ + Q_diff)]

    ptide = reduced_periods[0]
    lgQ_samples = get_dissipation_samples(prior_samples, tidal_period=10**ptide)

    q = scipy.stats.norm.cdf(1  )
    rl_stats = find_emcee_quantiles(lgQ_samples,q,0.001,1000)

    qunatile_val = rl_stats[0]

    

    lgQ_samples = lgQ_samples[max_burnins:].flatten()


    x,y = FFTKDE(kernel='gaussian', bw='silverman').fit(lgQ_samples).evaluate()

    b1=FFTKDE(kernel='gaussian', bw='silverman').fit(lgQ_samples).bw
    print('bw_silverman = {}'.format(b1))

    b2=FFTKDE(kernel='gaussian', bw='ISJ').fit(lgQ_samples).bw
    print('ISJ = {}'.format(b2))



    _dist = DiscreteSampling(lgQ_samples, b1)

    print(f'{ptide}\t{qunatile_val}\t{_dist.ppf(q)}')


if __name__=='__main__':

    args = cmdline_args()

    system_filename = args.chains_file
    print('\nTesting for {}'.format(system_filename))
    if system_filename is not None:
        system_kic = system_filename.split('.')[0].split('_')[1]

        prior_samples = corrected_samples(system_filename)
        
    # SAVE CONVERGENCE TEST FOR PERIOD [0.5,...,50]
    tidal_period_dependence(system_kic, prior_samples)
    
#-------------------------------------------------------------------------------
    #random codes might be useful at some point
#--------------------------------------------------------------------------------
    # Basically print the systems which still hasnt converged
    # all_systems = glob.glob('*.h5')
    # remaining_systems = []
    # for system_filename in all_systems:
    #     system_kic = system_filename.split('.')[0].split('_')[1]
    #     if system_kic not in converged_systems:
    #         remaining_systems.append(system_kic)
    # print(len(remaining_systems))
    # print(' '.join(remaining_systems))
#---------------------------------------------------------------------------------
    # Calculate convergence parameters at tidal period = 10days
    # save_filename = args.convergence_filename
    # if save_filename is not None:
    #     with open(save_filename,'a') as f:
    #         f.write(system_kic+'\n')
    # txt_file = open('period_dependence/new_stop_systems.txt')
    # stopped_systems = txt_file.read().split('\n')
    # if '' in stopped_systems:
    #     stopped_systems.remove('')
    # if system_kic not in stopped_systems:
    #     lgQ_samples = get_dissipation_samples(prior_samples)
    #     print(emcee_convergence_test(lgQ_samples, save_filename=save_filename))
    #     print('\n\n\n')
#---------------------------------------------------------------------------------
