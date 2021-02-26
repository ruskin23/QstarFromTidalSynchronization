import scipy
from scipy import integrate
from scipy import interpolate
from scipy.stats import norm
from scipy.misc import derivative
from utils import _get_filename,_get_chain,_fill_parameters,_cummulative_distribution
import numpy
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys
import pickle

def logQ_subplot(s,D,M):

    x=numpy.linspace(5,12,100000)

    # y=D['12']/max(D['12'])
    # plt.plot(x,M)
    # plt.plot(x,y)
    # plt.show()

    M_Z=interpolate.InterpolatedUnivariateSpline(x,M).integral(5,12)
    M=M/M_Z
    M=M/max(M)
    # s=numpy.reshape(numpy.array(s),(7,6))

    # _,axs=plt.subplots(7,6,sharex='col',sharey='row')
    plt.figure(figsize=(15,10))
    gs1=gridspec.GridSpec(6,7)
    gs1.update(wspace=0.0, hspace=0.0)

    for i in range(42):

        f_Z=interpolate.InterpolatedUnivariateSpline(x,D[s[i]]).integral(5,12)
        f=D[s[i]]/f_Z
        f=f/max(f)

        ax1=plt.subplot(gs1[i])
        plt.axis('on')
        ax1.plot(x,f)
        ax1.plot(x,M)
        # ax1.set_xticklabels([])
        # ax1.set_yticklabels([])
        ax1.tick_params(axis='x', labelsize= 5)
        ax1.tick_params(axis='y', labelsize= 5)
        ax1.label_outer()
        # ax1.set_aspect('equal')


    # for i in range(7):
    #     for j in range(6):
    #         f_Z=interpolate.InterpolatedUnivariateSpline(x,D[s[i][j]]).integral(5,12)
    #         f=D[s[i][j]]/f_Z
    #         f=f/max(f)
    #         print(s[i][j])
    #         print(len(x))
    #         print(len(f))
    #         axs[i,j].plot(x,f)
    #         axs[i,j].plot(x,M)
    #         axs[i,j].tick_params(axis='x', labelsize= 3)
    #         axs[i,j].tick_params(axis='y', labelsize= 3)

    # for ax in axs.flat:
    #     ax.label_outer()

    # plt.show()
    plt.savefig('all_pdf.png')

def logQ_M_plots(s,D,M):

    x=numpy.linspace(5,12,100000)

    M_Z=interpolate.InterpolatedUnivariateSpline(x,M).integral(5,12)
    M=M/M_Z
    M=M/max(M)

    

    for system in s:
        print(system)
        f_Z=interpolate.InterpolatedUnivariateSpline(x,D[system]).integral(5,12)
        f=D[system]/f_Z
        f=f/max(f)
        plt.plot(x,f)
        plt.plot(x,M)
        plt.savefig(f'pdf/compare/system_{system}.png')
        plt.close()





def get_samples(system,parameter):

    def get_chain():
        CHAIN=numpy.zeros(1)

        for c in ['ganymede','stampede']:
            for i in range(1,6):
                filename=_get_filename(system,c,i)
                filled_file=_fill_parameters(filename)
                chain=_get_chain(parameter,filled_file)
                if filename!=f'ganymede/MCMC_{system}/accepted_parameters_1.txt':
                    start_point=chain[0]
                    count_repeat=numpy.count_nonzero(chain==start_point)
                    chain=chain[count_repeat:]
                CHAIN=numpy.concatenate((CHAIN,chain),axis=None)
        CHAIN=CHAIN[1:]
        return CHAIN
    
    CHAIN=get_chain()
    cdf=_cummulative_distribution(CHAIN)
    Z=list(zip(*cdf))
    param_values=numpy.array(list(Z[0]))
    param_cdf=numpy.array(list(Z[1]))
    param_interp=interpolate.interp1d(param_values,param_cdf)
    param_reverse_interp=interpolate.interp1d(param_cdf,param_values)

    return param_values,param_cdf,param_interp,param_reverse_interp

def logQ_parameter_plots(s):
    parameters=['Porb','eccentricity','Wdisk','primary_mass','age','feh']

    for parameter in parameters:
        print(f'\n{parameter}')

        mass_error=[]
        logQ_error=[]

        mass_mean=[]
        logQ_mean=[]

        for system in s:        
            print(system)
            _,_,_,mass_reverse_interp=get_samples(system,parameter)
            _,_,_,logQ_reverse_interp=get_samples(system,'logQ')

            mass_lower_limit=mass_reverse_interp(0.3125)
            logQ_lower_limit=logQ_reverse_interp(0.3125)

            mass_upper_limit=mass_reverse_interp(0.6875)
            logQ_upper_limit=logQ_reverse_interp(0.6875)

            mass_error.append(mass_upper_limit-mass_lower_limit)
            logQ_error.append(logQ_upper_limit-logQ_lower_limit)

            mass_mean.append(mass_reverse_interp(0.5))
            logQ_mean.append(logQ_reverse_interp(0.5))

        plt.scatter(logQ_mean,mass_mean)
        plt.errorbar(logQ_mean,mass_mean,xerr=logQ_error,yerr=mass_error,linestyle='None')
        plt.xlabel('logQ')
        plt.ylabel(parameter)
        plt.savefig(f'pdf/logQ_{parameter}.png')
        plt.close()   

s=['1', '106', '109', '12', '120', '123', '126', '13', '137', '17', '20', '25', '28', '32', '36', '39', '43', '44', '47', '48', '50', '54', '56', '57', '67', '70', '73', '76', '79', '8', '80', '81', '83', '84', '85', '86', '88', '92', '93', '94', '95', '96']
with open('all_pdf_data.pickle','rb') as f:
    D=pickle.load(f)

M=D['M']
# logQ_M_plots(s,D,M)
logQ_subplot(s,D,M)

