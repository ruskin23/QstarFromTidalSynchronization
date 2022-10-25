import pickle
import sys
from tkinter.messagebox import NO
import matplotlib.pyplot as plt
import numpy

from pathlib import Path
from directories import directories

home_dir=str(Path.home())
path=directories(home_dir)
sys.path.append(path.poet_path+'/PythonPackage')
sys.path.append(path.poet_path+'/scripts')
sys.path.append('/home/ruskin/projects/')
sys.path.append('/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/version2_emcee/h5analysis/')
sys.path.append('/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/version2_emcee/')

from orbital_evolution.transformations import phase_lag, lgQ
from general_purpose_python_modules.emcee_quantile_convergence import *
from utils import *
from h5_analysis import *
import scipy
from KDEpy import FFTKDE
import json
from scipy import  interpolate


tex_fonts = {
    # Use LaTeX to write all text
    "text.usetex": True,
    "font.family": "serif",
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 14,
    "font.size": 10,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 8,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12
}

plt.rcParams.update(tex_fonts)


_quantiles=[scipy.stats.norm.cdf(c) for c in [-2,-1,1,2]]

def _normalization_constant(x, y, a, b):
    return interpolate.InterpolatedUnivariateSpline(x, y).integral(a, b)

def _kernel_bandwidth(samples):

    return FFTKDE(kernel='gaussian', bw='silverman').fit(samples).bw


class plot_distributions:


    def get_file_params(self):

        burnins = []
        
        with open(self.period_quantile_file,'r') as f:
            next(f)
            for idx, lines in enumerate(f):
                x = lines.split()
                #For burn-in
                for i in [3,6,9,12]:
                    burnins.append(int(x[i].split('/')[0]))

                #save params in file as grid
                self.QuantilePeriodGrid[idx] = numpy.array([[float(x[k]) for k in [0,1,4,7,10]]])

        self.burn_in = max(burnins)

        self.QuantilePeriodGrid = numpy.transpose(self.QuantilePeriodGrid)
        lgQ_2sigma_upper = self.QuantilePeriodGrid[-1]
        lgQ_2sigma_lower = self.QuantilePeriodGrid[1]
        self.MinDiffIdx = numpy.argmin(lgQ_2sigma_upper - lgQ_2sigma_lower)
        self.minQ = lgQ_2sigma_upper[self.MinDiffIdx]
        
        self.period_BestRange = self.QuantilePeriodGrid[0][lgQ_2sigma_upper < (self.minQ + self.Q_diff)]
        self.period_LeftBound, self.period_RightBound = self.period_BestRange[0], self.period_BestRange[-1]
        
        spin_filename = '/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/version2_emcee/catalog/filtering/nominal_value_catalog_Iconv_cutoff.txt'
        with open(spin_filename,'r') as f:
            next(f)
            for lines in f:
                x = lines.split()
                if x[1] == self.system_kic:
                    self.spin_dict['porb'] = float(x[2])
                    self.spin_dict['pspin'] = float(x[3])
                    self.spin_dict['error'] = float(x[4])


    def __init__(self,
                system_kic):

        self.system_kic = system_kic
        self.prior_samples = corrected_samples(f'system_{self.system_kic}.h5')

        self.period_array = numpy.linspace(numpy.log10(0.5), numpy.log10(50), 50)
        self.lgQArray = numpy.linspace(5, 20, 250)
        self.QuantilePeriodGrid = numpy.empty((50,5),float)

        self.period_quantile_file = f'period_dependence/{self.system_kic}.txt'
        self.burn_in = None
        self.period_LeftBound, self.period_RightBound = None, None
        self.period_BestRange = None
        self.Q_diff = 0.5
        self.minQ = None
        self.MinDiffIdx = None
        self.spin_dict = {}


        self.get_file_params()
    
    def get_pdfs(self):
        
        pdf_grid = numpy.empty((50,250), dtype=float)
        for idx, ptide in enumerate(self.period_array):
            if ptide < self.period_LeftBound and ptide < self.period_RightBound:
                pdf_grid[idx] = numpy.ones(250)/250.0
            else:
                lgQ_samples = get_dissipation_samples(self.prior_samples, tidal_period=ptide)[self.burn_in:].flatten()
                bw = _kernel_bandwidth(lgQ_samples)
                lgQ_dist = DiscreteSampling(lgQ_samples, bw)
                pdf_Q = lgQ_dist._pdf(self.lgQArray)
                pdf_Q /= _normalization_constant(self.lgQArray, pdf_Q, 5, 20)
                pdf_grid[idx] = pdf_Q
        return pdf_grid


    def plot_quantile(self, fig, ax, xlabel=False, ylabel=False):

        periods = self.QuantilePeriodGrid[0]
        ax.plot(periods, self.QuantilePeriodGrid[1], color='g', linewidth=1)
        ax.plot(periods, self.QuantilePeriodGrid[-1], color='g', linewidth=1)

        ax.plot(periods, self.QuantilePeriodGrid[2], color='k', linewidth=3)
        ax.plot(periods, self.QuantilePeriodGrid[3], color='k', linewidth=3)

        ax.set_title(r'KIC{} $Porb = {}$ $Pstar = {}\pm{}$'.format(self.system_kic, self.spin_dict['porb'], self.spin_dict['pspin'], self.spin_dict['error']))
        if xlabel: ax.set_xlabel(r'$\log_{10}{P_{tide}}$')
        if ylabel: ax.set_ylabel(r'$\log_{10}{Q_{\ast}^{\prime}}$')
        ax.set_ylim(5,16)
        # ax.set_xscale('log')
        # fig.savefig(f'test.png')
    
    def shade_region(self, fig, ax):
        
        ax.vlines(self.period_LeftBound, 5, 16, colors='k', linestyles='--' )
        ax.vlines(self.period_RightBound, 5, 16, colors='k', linestyles='--' )
        ax.fill_between(self.period_BestRange, self.QuantilePeriodGrid[-1][self.QuantilePeriodGrid[-1] < (self.minQ + self.Q_diff)], self.QuantilePeriodGrid[1][self.QuantilePeriodGrid[-1] < (self.minQ + self.Q_diff)])

def figure_1():

    with open('period_dependence/new_stop_systems.txt','r') as f:
        system_kic = f.read().split('\n')
    if '' in system_kic:
        system_kic.remove('')
    print(system_kic)
    # system_chunks = [system_kic[i:i+4] for i in range(0, len(system_kic), 4)]
    # print(system_chunks)

    # for k,chunk in enumerate(system_chunks):

    # fig = plt.figure(figsize=(15,15))
    # fig.subplots_adjust(hspace=0.5)
    # print(f'\nchunk = {k}')

    
    for i, kic in enumerate(system_kic):
        fig, ax = plt.subplots(1,1)
        print(f'system = {kic}')
        q = plot_distributions(kic)
        # ax = plt.subplots(1,1)
        q.plot_quantile(fig, ax, xlabel = True, ylabel = True)
        q.shade_region(fig, ax)
    # plt.savefig('/home/ruskin/projects/PhDDissertation2022/individual_constraints.pdf', bbox_inches='tight')
        plt.savefig(f'plots/individual_constraint_{kic}.pdf')
        plt.close()

def figure_2():

    with open('period_dependence/new_stop_systems.txt','r') as f:
        system_kic = f.read().split('\n')
    if '' in system_kic:
        system_kic.remove('')
    print(system_kic)
    with open('distributions_dict.json','r') as f:
        dist_dict = json.load(f)

    #Combined Distribution Grid Calculation
    individual_quantile_dict = {}
    pdf_combined_distribution = numpy.ones((50,250), dtype=float)

    for kic in system_kic:
        print(kic)
        q = plot_distributions(kic)
        pdf_grid = q.get_pdfs()
        
        #Adjust the pdf to whether to discard a limit or system
        #if discard == both, then the system will have uniform distribution
        bound_condition = dist_dict[kic]['discard']
        lgQ_best_pdf = pdf_grid[(q.period_array > q.period_LeftBound) & (q.period_array < q.period_RightBound)]
        print(bound_condition)
        for idx in range(len(lgQ_best_pdf)):
            pdf_lgQ = lgQ_best_pdf[idx]
            if bound_condition != 'None':
                if bound_condition == 'lower':
                    pdf_lgQ[:numpy.argmax(pdf_lgQ)] = max(pdf_lgQ)
                if bound_condition == 'upper':
                    pdf_lgQ[numpy.argmax(pdf_lgQ):] = max(pdf_lgQ)
                if bound_condition == 'both':
                    pdf_lgQ[:] = max(pdf_lgQ)
                pdf_lgQ /= _normalization_constant(q.lgQArray, pdf_lgQ, 5, 20)
            lgQ_best_pdf[idx] = pdf_lgQ
        pdf_grid[(q.period_array > q.period_LeftBound) & (q.period_array < q.period_RightBound)] = lgQ_best_pdf
        pdf_combined_distribution *= pdf_grid

        kic_dict = {}
        kic_dict['best_period'] = q.period_array[q.MinDiffIdx]
        lgQ_min_pdf = pdf_grid[q.MinDiffIdx]
        dist_object = DiscreteDistributions(lgQ_min_pdf, q.lgQArray, (5, 20))
        with open(q.period_quantile_file, 'r') as f:
            for lines in f:
                x = lines.split()
                if x[0] == str(q.period_array[q.MinDiffIdx]):
                    print(x)
                    break
        print(q.period_array[q.MinDiffIdx])
        print(q.burn_in)
        kic_dict['lowerlgQ_2'] = dist_object._ppf(scipy.stats.norm.cdf(-2))
        kic_dict['lowerlgQ_1'] = dist_object._ppf(scipy.stats.norm.cdf(-1))
        kic_dict['upperlgQ_1'] = dist_object._ppf(scipy.stats.norm.cdf(1))
        kic_dict['upperlgQ_2'] = dist_object._ppf(scipy.stats.norm.cdf(2))
        print(kic_dict)
        individual_quantile_dict[kic] = kic_dict
    
    #Plotting Figure Showing Combined Constraint
    print('plotting')
    fig, ax = plt.subplots(1,1)
    periods = numpy.linspace(numpy.log10(0.5), numpy.log10(50), 50)
    lgQ_array = numpy.linspace(5, 20, 250)
    combined_percentile_values = []
    
    for idx, ptide in enumerate(periods):
        lgQ_pdf = pdf_combined_distribution[idx]
        lgQ_pdf /= interpolate.InterpolatedUnivariateSpline(lgQ_array, lgQ_pdf).integral(5, 20)
        dist_object = DiscreteDistributions(lgQ_pdf, lgQ_array, (5, 20))
        combined_percentile_values.append(numpy.array([dist_object._ppf(p) for p in _quantiles]))
        
    for _, items_dict in individual_quantile_dict.items():
        ax.vlines(items_dict['best_period'], items_dict['lowerlgQ_1'], items_dict['upperlgQ_1'], linewidth = 3)
        ax.vlines(items_dict['best_period'], items_dict['lowerlgQ_2'], items_dict['upperlgQ_2'], linewidth = 1)
    
    combined_percentile_values = numpy.transpose(numpy.array(combined_percentile_values))
    upper_two_sigma = combined_percentile_values[3]
    for cpv in combined_percentile_values:
        ax.plot(periods[upper_two_sigma < 17], cpv[upper_two_sigma < 17], color = 'k')

    ax.set_xlabel(r'$\log_{10}{P_{tide}}$')
    ax.set_ylabel(r'$\log_{10}{Q_{\ast}^{\prime}}$')
    ax.set_yticks([k for k in range(5, 20, 2)])
    ax.set_xlim(numpy.log10(0.5), numpy.log10(50))

    with open('combined.npy', 'wb') as f:
        numpy.save(f, pdf_combined_distribution)
        numpy.save(f, combined_percentile_values)

    fig.savefig('plots/combined_constraints.pdf')
    # plt.savefig('/home/ruskin/projects/PhDDissertation2022/combined_constraints.pdf', bbox_inches='tight')
    plt.close()

    #Calculate Minimum Constraint
    fig2, ax2 = plt.subplots(1,1)
    min_idx = numpy.argmin(combined_percentile_values[-1] - combined_percentile_values[1])
    bestQ_pdf = pdf_combined_distribution[min_idx]
    dist_object = DiscreteDistributions(bestQ_pdf, lgQ_array, (5, 20))
    percentiles_bestQ = numpy.array([dist_object._ppf(p) for p in _quantiles])
    median_Q = dist_object._ppf(0.5)
    bestQ_pdf /= max(bestQ_pdf)
    ax2.plot(lgQ_array, bestQ_pdf, color = 'k', label=r'Frequency-Dependent $Q_{\ast}^{\prime}$')


    lgQ_array_old = numpy.linspace(5, 12, 100000)
    with open('/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/version1_metropolis_hasting/SolutionFileBreaks0.0.txt', 'r') as f:
        next(f)
        common_pdf = 1
        for lines in f:
            x = lines.split()
            system_number = x[0]
            with open('all_pdf_data.pickle','rb') as f:
                D=pickle.load(f)
            if system_number not in ['31', '57']:
                f_Z = interpolate.InterpolatedUnivariateSpline(lgQ_array_old, D[system_number]).integral(5,12)
                f = D[system_number]/f_Z
                common_pdf *= f
        common_pdf /= interpolate.InterpolatedUnivariateSpline(lgQ_array_old, common_pdf).integral(5,12)
        common_pdf /= max(common_pdf)

    ax2.plot(lgQ_array_old, common_pdf, color = 'r', label = r'Constant $Q_{\ast}^{\prime}$')
    ax2.legend()
    ax2.set_xlim((5,12))
    

    plt.savefig('plots/comparison.pdf', bbox_inches='tight')
    plt.close()
    print(percentiles_bestQ, median_Q)


def figure_3():

    with open('combined.npy', 'rb') as f:
        pdf_combined_distribution = numpy.load(f)
        combined_percentile_values = numpy.load(f)


    #Calculate Minimum Constraint
    fig2 = plt.figure(figsize=(5, 5))
    ax2 = fig2.add_subplot(111)
    lgQ_array = numpy.linspace(5, 20, 250)
    min_idx = numpy.argmin(combined_percentile_values[-1] - combined_percentile_values[1])
    bestQ_pdf = pdf_combined_distribution[min_idx]
    dist_object = DiscreteDistributions(bestQ_pdf, lgQ_array, (5, 20))
    percentiles_bestQ = numpy.array([dist_object._ppf(p) for p in _quantiles])
    median_Q = dist_object._ppf(0.5)
    bestQ_pdf /= max(bestQ_pdf)
    ax2.plot(lgQ_array, bestQ_pdf, color = 'k', label=r'Frequency-Dependent $Q_{\ast}^{\prime}$')


    lgQ_array_old = numpy.linspace(5, 12, 100000)
    with open('/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/version1_metropolis_hasting/SolutionFileBreaks0.0.txt', 'r') as f:
        next(f)
        common_pdf = 1
        for lines in f:
            x = lines.split()
            system_number = x[0]
            with open('all_pdf_data.pickle','rb') as f:
                D=pickle.load(f)
            if system_number not in ['31', '57']:
                f_Z = interpolate.InterpolatedUnivariateSpline(lgQ_array_old, D[system_number]).integral(5,12)
                f = D[system_number]/f_Z
                common_pdf *= f
        common_pdf /= interpolate.InterpolatedUnivariateSpline(lgQ_array_old, common_pdf).integral(5,12)
        common_pdf /= max(common_pdf)

    ax2.plot(lgQ_array_old, common_pdf, color = 'r', label = r'Constant $Q_{\ast}^{\prime}$')
    ax2.legend()
    ax2.set_xlim((5,12))
    ax2.set_xlabel(r"$\log_{10}{}Q_{\ast}^{\prime}$")
    ax2.set_ylabel('pdf')

    

    plt.savefig('/home/ruskin/projects/PhDDissertation2022/comparison_common.pdf', bbox_inches='tight')
    plt.close()
    print(percentiles_bestQ, median_Q)




    # with open('period_dependence/new_stop_systems.txt','r') as f:
    #     system_kic = f.read().split('\n')
    # if '' in system_kic:
    #     system_kic.remove('9971475')
    #     system_kic.remove('')
    
    # old_system_numbers = []
    # with open('/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/version1_metropolis_hasting/SolutionFileBreaks0.0.txt', 'r') as f:
    #     next(f)
    #     for lines in f:
    #         old_system_numbers.append(lines.split()[0])
    # common_kic = []
    # common_numbers = []
    # with open('/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/version1_metropolis_hasting/SpinlogQCatalog_el0.4.txt', 'r') as f:
    #     next(f)
    #     for lines in f:
    #         x = lines.split()
    #         number = x[0]
    #         kic = x[1]
    #         if number in old_system_numbers:
    #             if kic in system_kic:
    #                 common_kic.append(kic)
    #                 common_numbers.append(number)
    # print(common_kic, common_numbers)
    # for i, kic in enumerate(common_kic):
    #     print(kic)
    #     fig, ax = plt.subplots(1,1)
    #     q = plot_distributions(kic)
    #     new_pdf = q.get_pdfs()[q.MinDiffIdx]
    #     new_pdf /= max(new_pdf)
    #     ax.plot(q.lgQArray, new_pdf, color='k', label=r'Frequency-Dependent $Q_{\ast}^{\prime}$')

    #     with open('all_pdf_data.pickle','rb') as f:
    #         D=pickle.load(f)
        
    #     x = numpy.linspace(5, 12, 100000)
    #     f_Z = interpolate.InterpolatedUnivariateSpline(x, D[common_numbers[i]]).integral(5,12)
    #     f = D[common_numbers[i]]/f_Z
    #     f /= max(f)
    #     ax.plot(x, f, color='r', label=r'Constant $Q_{\ast}^{\prime}$')

    #     ax.set_xlim((5,12))
    #     ax.set_xlabel(r"$\log_{10}{}Q_{\ast}^{\prime}$")
    #     ax.set_ylabel('pdf')

    #     ax.set_title(f'KIC{kic}')
    #     ax.legend()

    #     fig.savefig(f'comparison_{kic}.png')
    #     plt.close()



if __name__ == '__main__':

    figure_1()
    figure_2()
    figure_3()
    # with open('combined.npy', 'rb') as f:
    #     pdf_combined_distribution = numpy.load(f)
    #     combined_percentile_values = numpy.load(f)
    
    # upper_two_sigma = combined_percentile_values[3]
    # print(upper_two_sigma)
    # periods = numpy.linspace(numpy.log10(0.5), numpy.log10(50), 50)
    # for cpv in combined_percentile_values:
    #     plt.plot(periods[upper_two_sigma < 17], cpv[upper_two_sigma < 17])
    # plt.savefig('test.png')

    # with open('all_pdf_data.pickle','rb') as f:
    #     D=pickle.load(f)
    
    # with open('/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/version1_metropolis_hasting/SpinlogQCatalog_el0.4.txt', 'r') as f:
    #     next(f)
    #     for lines in f:
    #         x = lines.split()
    #         if kic == x[1]:
    #             system_number = x[0]
    # x = numpy.linspace(5, 12, 100000)
    # f_Z = interpolate.InterpolatedUnivariateSpline(x, D[system_number]).integral(5,12)
    # f = D[system_number]/f_Z
    # f /= max(f)

    
    #dump dictionaries to json. 
    # data_dict = {}
    # with open('/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/version2_emcee/catalog/filtering/nominal_value_catalog_Iconv_cutoff.txt') as f:
    #     next(f)
    #     for line in f:
    #         x=line.split()
            
    #         kic_dict = {}
    #         kic_dict['orbital_period'] = x[2]
    #         kic_dict['spin_period'] = x[3]
    #         kic_dict['spin_error'] = x[4]
    #         kic_dict['eccentricity'] = x[5]
    #         kic_dict['discard'] = 'None'

    #         data_dict[x[1]] = kic_dict

    # with open('distributions_dict.json','w') as f:
    #     json.dump(data_dict, f, indent=4)