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

from orbital_evolution.transformations import phase_lag, lgQ
from general_purpose_python_modules.emcee_quantile_convergence import *
import scipy

_quantities=['primary_mass',
             'secondary_mass',
             'feh',
             'age',
             'eccentricity',
             'phase_lag_max',
             'alpha',
             'omegeref']


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

def cmdline_args():

    p = ArgumentParser()
    p.add_argument('-s',
                    dest = 'chains_file')

    p.add_argument('-p',
                   dest = 'parameter')

    p.add_argument('-c',
                    dest='convergence_filename',
                    default=None)

    p.add_argument('--all',
                    action='store_true',
                    default=False)

    p.add_argument('--steps',
                   action='store_true',
                   default=False)

    return p.parse_args()



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


def emcee_convergence_test(lgQ_samples,
                           save_filename=None):

    total_steps = numpy.shape(lgQ_samples)[0]
    _quantities=[scipy.stats.norm.cdf(c) for c in [-2,-1,1,2]]
    st_list = []
    for qunatile in _quantities:
        print('Calculating for quantile {}'.format(qunatile))
        a = find_emcee_quantiles(lgQ_samples,qunatile,0.001,1000)
    
        qunatile_val = a[0]
        r = a[2]
        thin = a[3]
        burnin = a[4]
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
            print(a)
            st_list.append(a)


    return st_list
    
def tidal_period_dependence(system_kic, prior_samples):

    tidal_period = numpy.linspace(0.5,50,100)
    quantiles = [scipy.stats.norm.cdf(c) for c in [-2,-1,1,2]]
    f = open('period_dependence/'+system_kic+'.txt','w')
    f.write('period'+'\t'+
            '-2sigma({})'.format(repr(quantiles[0]))+'\t'+
            'r'+'\t'+
            'burnin/total_steps'+'\t'+
            '-1sigma({})'.format(repr(quantiles[0]))+'\t'+
            'r'+'\t'+
            'burnin/total_steps'+'\t'+
            '1sigma({})'.format(repr(quantiles[0]))+'\t'+
            'r'+'\t'+
            'burnin/total_steps'+'\t'+
            '2sigma({})'.format(repr(quantiles[0]))+'\t'+
            'r'+'\t'+
            'burnin/total_steps'+'\n')

    for p in tidal_period:
        lgQ_samples = get_dissipation_samples(prior_samples, tidal_period=p)
        st = emcee_convergence_test(lgQ_samples)

        l = '\t'.join(['\t'.join([repr(s[0]), repr(s[2]), '{}/{}'.format(repr(s[4]),repr(len(lgQ_samples)))]) for s in st])

        f.write(repr(p)+'\t'+l+'\n')
    f.close()

def get_dissipation_samples(samples, tidal_period=10):

    q_samples = []
    for s in samples:
        delta0 = s[5]
        alpha = s[6]
        omega_br = s[7]
        omega_min = 2*numpy.pi/50
        omega = 2*numpy.pi/tidal_period
        
        if alpha > 0:
            if omega < omega_min:
                delta = (omega_min/omega_br)**alpha
            elif omega > omega_min and omega < omega_br:
                delta = (omega/omega_br)**alpha
            else:
                delta = 1
        else:
            if omega < omega_br:
                delta = 1
            else:
                delta = (omega/omega_br)**alpha
        
        q_samples.append(lgQ(delta0*delta))

    if numpy.nan in q_samples:        
        print('NaN found in Q_SAMPLES!!!!!!!!!')
    q_samples = numpy.reshape(q_samples,(len(q_samples)//64,64))

    return q_samples


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





if __name__=='__main__':

    # data = numpy.empty((100,5), float)
    # with open('period_dependence/4678171.txt','r') as f:
    #     next(f)
    #     for i,lines in enumerate(f):
    #         x=lines.split()
    #         data[i] = numpy.array([numpy.log10(float(x[k])) if k == 0 else float(x[k]) for k in [0,1,4,7,10]])
    
    # data = numpy.transpose(data)
    # for i in range(4):
    #     plt.plot(data[0],data[i+1])
    # plt.savefig('test.png')

    args = cmdline_args()

    system_filename = args.chains_file

    save_filename = args.convergence_filename
    system_kic = system_filename.split('.')[0].split('_')[1]
    if save_filename is not None:
        
        with open(save_filename,'a') as f:
            f.write(system_kic+'\n')

    print('\nTesting for {}'.format(system_filename))


    prior_samples = corrected_samples(system_filename)
    # tidal_period_dependence(system_kic, prior_samples)
    lgQ_samples = get_dissipation_samples(prior_samples)
    print(emcee_convergence_test(lgQ_samples, save_filename=save_filename))

    print('\n\n\n')



