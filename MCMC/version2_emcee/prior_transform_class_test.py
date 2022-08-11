# from os import plock
import numpy
import time
from utils import erf_fun
import random
import corner
import matplotlib.pyplot as plt
import logging
import pickle

from scipy.interpolate import interp1d

from pathlib import Path
from directories import directories

home_dir=str(Path.home())
path=directories(home_dir)


# _logger = logging.getLogger(__name__)



class prior_transform:

    def __init__(self,system_num,bandwidth):

        self.system_num=system_num

        system_chains=numpy.load(path.current_directory+'/catalog/samples/chains/'+system_num+'.npz')
        all_chains=numpy.transpose(system_chains['thinned_chain'])

        self.M_samples=all_chains[0]
        self.Q_samples=all_chains[1]
        self.Z_samples=all_chains[2]
        self.t_samples=all_chains[3]

        ecosw_samples=all_chains[8]
        esinw_samples=all_chains[9]
        self.e_samples=numpy.sqrt(ecosw_samples**2 + esinw_samples**2)

        minmaxdiff_Z=abs(max(self.Z_samples)-min(self.Z_samples))
        minmaxdiff_e=abs(max(self.e_samples)-min(self.e_samples))
        self.M_prior=numpy.linspace(min(self.M_samples), max(self.M_samples), 10000)
        self.Q_prior=numpy.linspace(min(self.Q_samples), max(self.Q_samples), 10000)
        self.Z_prior=numpy.linspace(min(self.Z_samples)-minmaxdiff_Z, max(self.Z_samples)+minmaxdiff_Z,10000)
        # self.Z_prior=numpy.linspace(min(self.Z_samples), max(self.Z_samples),10000)
        self.t_prior=numpy.linspace(min(self.t_samples), max(self.t_samples), 10000)
        self.e_prior=numpy.linspace(max(min(self.e_samples),0),min(max(self.e_samples),1),10000)
        print(self.e_prior)

        self.bandwidth=bandwidth

        self.h_M=self.bandwidth[0]
        self.h_Q=self.bandwidth[1]
        self.h_Z=self.bandwidth[2]
        self.h_t=self.bandwidth[3]
        self.h_e=self.bandwidth[4]
        # self.h_e=0.0002


    def paramter_evaluate(self,uniform_values,Z_cdf_interp,erfq,erfm,erft,erfe):
    # def paramter_evaluate(self,uniform_values):

        # Z_cdf=0.0
        # for s in self.Z_samples:
        #     z=(self.Z_prior-s)/self.h_Z
        #     Z_cdf+=erf_fun(z)
        # Z_cdf/=len(self.Z_samples)
        # Z_cdf_interp = interp1d(Z_cdf,self.Z_prior)

        # q=(self.Q_prior[None,:]-self.Q_samples[:,None])/self.h_Q
        # m=(self.M_prior[None,:]-self.M_samples[:,None])/self.h_M
        # t=(self.t_prior[None,:]-self.t_samples[:,None])/self.h_t

        # erfq=erf_fun(q)
        # erfm=erf_fun(m)
        # erft=erf_fun(t)

        # if uniform_values[0]<min(Z_cdf):Z_value=min(self.Z_samples)
        # elif uniform_values[0]>max(Z_cdf):Z_value=max(self.Z_samples)
        # else:Z_value=Z_cdf_interp(uniform_values[0])
        Z_value=Z_cdf_interp(uniform_values[0])
        Z_p=numpy.exp(-0.5*(((Z_value-self.Z_samples)/self.h_Z)**2))

        F_Q=numpy.dot(Z_p,erfq)
        F_Q/=F_Q.max()
        Q_value=self.Q_prior[numpy.argmin(abs(F_Q-uniform_values[1]))]
        Q_p=numpy.exp(-0.5*(((Q_value-self.Q_samples)/self.h_Q)**2))

        F_M=numpy.dot(Z_p*Q_p,erfm)
        F_M/=F_M.max()
        M_value=self.M_prior[numpy.argmin(abs(F_M-uniform_values[2]))]
        M_p=numpy.exp(-0.5*(((M_value-self.M_samples)/self.h_M)**2))

        F_t=numpy.dot(Z_p*Q_p*M_p,erft)
        F_t/=F_t.max()
        t_value=self.t_prior[numpy.argmin(abs(F_t-uniform_values[3]))]
        t_p=numpy.exp(-0.5*(((t_value-self.t_samples)/self.h_t)**2))


        F_e=numpy.dot(Z_p*Q_p*M_p*t_p,erfe)
        F_e/=F_e.max()
        e_value=self.t_prior[numpy.argmin(abs(F_e-uniform_values[3]))]


        # e_value=e_cdf_interp(uniform_values[4])
        # count, bins_count = numpy.histogram(self.e_samples, bins=100)
        # pdf=count/sum(count)
        # e_cdf=numpy.cumsum(pdf)
        # e_cdf_interp = interp1d(e_cdf,self.e_prior)
        # if uniform_values[4]<min(e_cdf):e_value=min(self.e_samples)
        # elif uniform_values[4]>max(e_cdf):e_value=max(self.e_samples)
        # else:e_value=e_cdf_interp(uniform_values[4])        
        # e_value=self.e_prior[numpy.argmin(abs(e_cdf-uniform_values[4]))]

        #---------------------------------------------------------------------------->
        return [Z_value,M_value,Q_value,t_value,e_value] 
        #<------------------------------------------------------------------------------


def corner_plot(system_number,bandwidth=None,method=None):

    prior=prior_transform(system_number,bandwidth)

    Z_cdf=0.0
    for s in prior.Z_samples:
        z=(prior.Z_prior-s)/prior.h_Z
        Z_cdf+=erf_fun(z)
    Z_cdf/=len(prior.Z_samples)
    Z_cdf_interp = interp1d(Z_cdf,prior.Z_prior)

    # e_cdf=0.0
    # for s in prior.e_samples:
    #     e=(prior.e_prior-s)/prior.h_e
    #     e_cdf+=erf_fun(e)
    # e_cdf/=len(prior.e_samples)
    # e_cdf_interp = interp1d(e_cdf,prior.e_prior)


    q=(prior.Q_prior[None,:]-prior.Q_samples[:,None])/prior.h_Q
    m=(prior.M_prior[None,:]-prior.M_samples[:,None])/prior.h_M
    t=(prior.t_prior[None,:]-prior.t_samples[:,None])/prior.h_t
    e=(prior.e_prior[None,:]-prior.e_samples[:,None])/prior.h_e

    erfq=erf_fun(q)
    erfm=erf_fun(m)
    erft=erf_fun(t)
    erfe=erf_fun(e)

    values=[]
    for i in range(5000):
        rnd_nums=[random.uniform(min(Z_cdf),max(Z_cdf)),random.uniform(0,1),random.uniform(0,1),random.uniform(0,1),random.uniform(0,1)]
        values.append(prior.paramter_evaluate(rnd_nums,Z_cdf_interp,erfq,erfm,erft,erfe))

    flat_values = [item for sublist in values for item in sublist]

    d=numpy.array(numpy.vstack(numpy.transpose(numpy.array(flat_values))))
    data_corner=d.reshape([len(d)//5,5])
    figure=corner.corner(data_corner)
    # plt.show()
    plt.savefig('compare_plots_bandwidths/system_{}/corner_{}_{}.png'.format(system_number,method,system_number))


def get_bandwidth(param,data,plot=True,get_value=True,method=None):

    print('\nParam: {}'.format(param))

    x,y = FFTKDE(kernel='gaussian', bw='silverman').fit(data).evaluate()
    y = FFTKDE(kernel='gaussian', bw='ISJ').fit(data).evaluate(x)

    if method == 'silverman':
        bwdth=FFTKDE(kernel='gaussian', bw='silverman').fit(data).bw
        print('bw_silverman = {}'.format(bwdth))
        # plt.plot(x, y, label='KDE /w silverman')

    if method == 'ISJ':
        bwdth=FFTKDE(kernel='gaussian', bw='ISJ').fit(data).bw
        print('ISJ = {}'.format(bwdth))
        # plt.plot(x, y, label='KDE /w ISJ')
    
    if plot==True:
        plt.grid(True, ls='--', zorder=-15); plt.legend()
        plt.show()
    if get_value:
        return bwdth


# '9892471'


def get_sampling_utils(system_number):

    with open(path.current_directory+'/bandwidth.pickle','rb') as f:
        bw=pickle.load(f)
    bandwidths=bw[system_number]['ISJ']
    h_M=bandwidths[0]
    h_Q=bandwidths[1]
    h_Z=bandwidths[2]
    h_t=bandwidths[3]
    with open(path.current_directory+'/e_bandwidth.pickle','rb') as f:
        h_e=pickle.load(f)
    h_e=h_e[system_number]
    
    system_chains=numpy.load(path.current_directory+'/catalog/samples/chains/'+system_number+'.npz')
    all_chains=numpy.transpose(system_chains['thinned_chain'])

    M_samples=all_chains[0]
    Q_samples=all_chains[1]
    Z_samples=all_chains[2]
    t_samples=all_chains[3]

    ecosw_samples=all_chains[8]
    esinw_samples=all_chains[9]
    e_samples=numpy.sqrt(ecosw_samples**2 + esinw_samples**2)

    M_prior=numpy.linspace(min(M_samples), max(M_samples), 10000)
    Q_prior=numpy.linspace(min(Q_samples), max(Q_samples), 10000)
    Z_prior=numpy.linspace(min(Z_samples), max(Z_samples),10000)
    t_prior=numpy.linspace(min(t_samples), max(t_samples), 10000)
    e_prior=numpy.linspace(min(e_samples),max(e_samples),10000)

    q=(Q_prior[None,:]-Q_samples[:,None])/h_Q
    m=(M_prior[None,:]-M_samples[:,None])/h_M
    t=(t_prior[None,:]-t_samples[:,None])/h_t
    e=(e_prior[None,:]-e_samples[:,None])/h_e

    erfq=erf_fun(q)
    erfm=erf_fun(m)
    erft=erf_fun(t)
    erfe=erf_fun(e)

    sampling_util_dict=dict()
    sampling_util_dict['M_samples']=M_samples
    sampling_util_dict['Q_samples']=Q_samples
    sampling_util_dict['Z_samples']=Z_samples
    sampling_util_dict['t_samples']=t_samples
    sampling_util_dict['e_samples']=e_samples

    sampling_util_dict['M_prior']=M_prior
    sampling_util_dict['Q_prior']=Q_prior
    sampling_util_dict['Z_prior']=Z_prior
    sampling_util_dict['t_prior']=t_prior
    sampling_util_dict['e_prior']=e_prior

    sampling_util_dict['erfq']=erfq
    sampling_util_dict['erfm']=erfm
    sampling_util_dict['erft']=erft
    sampling_util_dict['erfe']=erfe

    sampling_util_dict['h_M']=h_M
    sampling_util_dict['h_Q']=h_Q
    sampling_util_dict['h_Z']=h_Z
    sampling_util_dict['h_t']=h_t
    sampling_util_dict['h_e']=h_e

    return sampling_util_dict

if __name__ == '__main__':

    UTILS=dict()
    with open(path.current_directory+'/catalog/filtering/nominal_value_catalog_Iconv_cutoff.txt','r') as f:
        next(f)
        for lines in f:
            x=lines.split()
            system_number=x[1]
            UTILS[system_number]=get_sampling_utils(system_number)

    with open(path.current_directory+'/SAMPLING_UTILS.pickle','wb') as f:
        pickle.dump(UTILS,f)

    # Bandwidths=dict()
    # with open('catalog/filtering/nominal_value_catalog_Iconv_cutoff.txt','r') as f:
    #     next(f)
    #     for lines in f:
    #         x=lines.split()
    #         system_number = '9892471'#x[1]
    #         print('\nSystem = {}'.format(system_number))
    #         system_chains=numpy.load(path.current_directory+'/catalog/samples/chains/'+system_number+'.npz')
    #         all_chains=numpy.transpose(system_chains['thinned_chain'])
    #         temp_dict=dict()
    #         for b in ['ISJ']:
    #             print('Bandwidth Method = {}'.format(b))
    #             bandwidth=[]
    #             params=['M','Q','Z','t','e']
    #             for i in range(4):
    #                 data=all_chains[i]
    #                 bandwidth.append(get_bandwidth(params[i],data,plot=False,get_value=True,method=b))
    #             ecosw_samples=all_chains[8]
    #             esinw_samples=all_chains[9]
    #             data=numpy.sqrt(ecosw_samples**2 + esinw_samples**2)
    #             bandwidth.append(get_bandwidth('e',data,plot=False,get_value=True,method=b))
    #             temp_dict[b]=bandwidth
    #             corner_plot(system_number,bandwidth=bandwidth,method=b)
    #         Bandwidths[system_number]=temp_dict
    #         break
    
    # with open('bandwidth2.pickle','wb') as f:
    #     pickle.dump(Bandwidths,f)

    
    # extra=False
    # if extra:
    #     for i in range(1):

    #         uniform_values=[random.uniform(0,1),random.uniform(0,1),random.uniform(0,1),random.uniform(0,1),random.uniform(0,1)]

    #         M_samples=all_chains[0]
    #         Q_samples=all_chains[1]
    #         Z_samples=all_chains[2]
    #         t_samples=all_chains[3]

    #         ecosw_samples=all_chains[8]
    #         esinw_samples=all_chains[9]
    #         e_samples=numpy.sqrt(ecosw_samples**2 + esinw_samples**2)

    #         M_prior=numpy.linspace(0.2, 2.5, 1000)
    #         Q_prior=numpy.linspace(0.0085, 1, 1000)
    #         Z_prior=numpy.linspace(0.001, 0.06, 1000)
    #         t_prior=numpy.linspace(6, 10.1, 1000)
    #         e_prior=numpy.linspace(0,1,1000)

    #         plot=False
    #         get_bandwidth('M',M_samples,M_prior,plot=plot)
    #         get_bandwidth('Q',Q_samples,Q_prior,plot=plot)
    #         get_bandwidth('Z',Z_samples,Z_prior,plot=plot)
    #         get_bandwidth('t',t_samples,t_prior,plot=plot)


    #         bandwidth=[0.01,0.001,0.0001,0.001]

    #         h_M=bandwidth[0]
    #         h_Q=bandwidth[1]
    #         h_Z=bandwidth[2]
    #         h_t=bandwidth[3]
    #         h_e=0.01

    #         timex=time.time()
    #         Z_cdf=0.0
    #         for s in Z_samples:
    #             z=(Z_prior-s)/h_Z
    #             Z_cdf+=erf_fun(z)
    #         Z_cdf/=len(Z_samples)
    #         Z_value=Z_prior[numpy.argmin(abs(Z_cdf-uniform_values[0]))]
    #         Z_p=numpy.exp(-0.5*(((Z_value-Z_samples)/h_Z)**2))
    #         print(time.time()-timex)
            
    #         q=(Q_prior[None,:]-Q_samples[:,None])/h_Q
    #         F_Q=numpy.dot(Z_p,q)
    #         F_Q/=F_Q.max()
    #         Q_value=Q_prior[numpy.argmin(abs(F_Q-uniform_values[1]))]
    #         Q_p=numpy.exp(-0.5*(((Q_value-Q_samples)/h_Q)**2))

    #         m=(M_prior[None,:]-M_samples[:,None])/h_M
    #         F_M=numpy.dot(Z_p*Q_p,m)
    #         F_M/=F_M.max()
    #         M_value=M_prior[numpy.argmin(abs(F_M-uniform_values[2]))]
    #         M_p=numpy.exp(-0.5*(((M_value-M_samples)/h_M)**2))

    #         t=(t_prior[None,:]-t_samples[:,None])/h_t
    #         F_t=numpy.dot(Z_p*Q_p*M_p,t)
    #         F_t/=F_t.max()
    #         t_value=t_prior[numpy.argmin(abs(F_t-uniform_values[3]))]

    # print(time.time()-time1)
    # system_chains=numpy.load(path.current_directory+'/catalog/samples/chains/'+system_number+'.npz')
    # all_chains=numpy.transpose(system_chains['thinned_chain'])

    # ecosw_samples=all_chains[8]
    # esinw_samples=all_chains[9]
    # e_samples=numpy.sqrt(ecosw_samples**2 + esinw_samples**2)
    # e_prior=numpy.linspace(0,1,1000)

    # h_e=0.01

    # e_cdf=0
    # for s in e_samples:
    #     alpha_e=(e_prior-s)/h_e
    #     e_cdf=e_cdf+erf_fun(alpha_e)
    # e_cdf=e_cdf/len(e_samples)
    # plt.plot(e_prior,e_cdf)
    # plt.xlim(0,0.1)

    # for bins in [100]:
    #     count, bins_count = numpy.histogram(e_samples, bins=bins)
    #     pdf=count/sum(count)
    #     cdf=numpy.cumsum(pdf)
    #     plt.plot(bins_count[1:],pdf)
    #     plt.plot(bins_count[1:],cdf)
    # plt.show()
    # rnd_nums=[random.uniform(0,1),random.uniform(0,1),random.uniform(0,1),random.uniform(0,1),random.uniform(0,1)]


    # sampler=prior_transform(system_number)

    # print(sampler.paramter_evaluate(rnd_nums))


    # grid = GridSearchCV(KernelDensity(),
    #                     {'bandwidth': numpy.linspace(0.1, 1.0, 30)},
    #                     cv=20) # 20-fold cross-validation
    # grid.fit(data[:, None])
    # print(grid.best_params_)

    # e_cdf=0
    # for s in self.e_samples:
    #     alpha_e=(self.e_prior-s)/self.h_e
    #     e_cdf=e_cdf+erf_fun(alpha_e)
    # e_cdf=e_cdf/len(self.e_samples)
