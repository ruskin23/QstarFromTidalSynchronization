import sys
import numpy 
import scipy
import pickle
import matplotlib.pyplot as plt
from scipy import integrate
from scipy import interpolate
from scipy.stats import norm
from scipy.optimize import curve_fit
from utils import _get_filename,_get_chain,_fill_parameters,_cummulative_distribution,adjust_chain



class p_value():

    def __init__(self,M_mean,M_interp,M_Z):

        
        self.clusters=['ganymede','stampede']

        self.logQ_array=numpy.linspace(5,12,100000)
        self.two_sided=False
        self.M_mean=M_mean
        self.M_interp=M_interp
        self.M_Z=M_Z
        

    def get_chain(self,system):

        # self.CHAIN=numpy.zeros(1)
        
        # for c in self.clusters:
        #     for i in range(1,6):
        #         filename=_get_filename(system,c,i)
        #         filled_file=_fill_parameters(filename)
        #         chain=_get_chain('logQ',filled_file)
        #         chain=adjust_chain(system,chain)
        #         if filename!=f'ganymede/MCMC_{system}/accepted_parameters_1.txt':
        #             start_point=chain[0]
        #             count_repeat=numpy.count_nonzero(chain==start_point)
        #             chain=chain[count_repeat:]
        #         self.CHAIN=numpy.concatenate((self.CHAIN,chain),axis=None)
        # self.CHAIN=adjust_chain(system,self.CHAIN[1:])

        with open('../complete_chains.pickle','rb') as f:
            D=pickle.load(f)
        for system_name,params in D.items():
            if system_name==system:
                for name,values in params.items():
                    if name=='logQ':self.CHAIN=values


    def get_samples(self,system):

        self.get_chain(system)

        cdf=_cummulative_distribution(self.CHAIN)
        Z=list(zip(*cdf))
        self.logQ_values=numpy.array(list(Z[0]))
        self.logQ_cdf=numpy.array(list(Z[1]))
        self.logQ_median=numpy.median(self.logQ_values)
        self.logQ_mean=numpy.mean(self.logQ_values)
        self.logQ_interp=interpolate.interp1d(self.logQ_values,self.logQ_cdf)
        self.logQ_reverse_interp=interpolate.interp1d(self.logQ_cdf,self.logQ_values)

    def get_cdf_statistic(self,quantile=None,value=None):
        if quantile is not None:
            return self.logQ_reverse_interp(quantile)
        
        if value is not None:
            if value<min(self.logQ_values):return 0
            elif value>max(self.logQ_values): return 1
            else:return self.logQ_interp(value)
    
    def logQ_pdf(self,x,h=None):

        N=len(self.CHAIN)
        if h is None:h=3.5*numpy.power(N,-1.0/3)

        f_h=0
        for i in range(N):
            y=(x-self.CHAIN[i])/h
            f_h=f_h+norm.pdf(y)

        return f_h/(N*h)

    def logQ_cdf_p(self,x):

        N=len(self.CHAIN)
        h=3.5*numpy.power(N,-1.0/3)

        f_h=0
        for i in range(N):
            y=(x-self.CHAIN[i])/h
            f_h=f_h+norm.cdf(y)

        return f_h/N

    def p_E(self):

        q_stat=self.get_cdf_statistic(value=self.M_mean)
        return 2*min(q_stat,1-q_stat)

    def integrand_E_p(self,x):

        cdf_p=self.get_cdf_statistic(value=x)
        q=2*min(cdf_p,1-cdf_p)
        return (1.0/self.M_Z)*q*self.M_interp(x)

    def E_p(self):

        p=integrate.quad(self.integrand_E_p,5,12)[0]
        return p
    
    def plot_check(self):

        I=numpy.zeros(len(self.logQ_array))
        for i,q in enumerate(self.logQ_array):
            I[i]=self.integrand_E_p(q)/(self.M_Z)

        I_interp=interpolate.InterpolatedUnivariateSpline(self.logQ_array,I)
        print(integrate.quad(I_interp,5,12))
        print(f'Spline Integration = {I_interp.integral(5,12)}')
        plt.plot(self.logQ_array,I)
        plt.show()

    def __call__(self,system,method):

        
        self.get_samples(system)
        if method=='p_E':return self.p_E()
        elif method=='E_p':return self.E_p()
        else: raise ValueError

def find_nearest(array, value):
    idx = (numpy.abs(array - value)).argmin()
    return idx

def gaussian_fit(x,mu,sigma):

    y=(x-mu)/sigma
    N=(1/(sigma*numpy.sqrt(2*numpy.pi)))

    return N*numpy.exp(-0.5*y*y)

def M_stats_pdf(logQ_array,M):

    mean=7.81759817598176
    sigma=0.10696106961069685
    value_l=mean-2*sigma
    value_r=mean+2*sigma
    logQ_idx_l=find_nearest(logQ_array,value_l)
    logQ_idx_r=find_nearest(logQ_array,value_r)

    M_new=M[logQ_idx_l:logQ_idx_r]
    logQ_new=logQ_array[logQ_idx_l:logQ_idx_r]

    plt.plot(logQ_new,M_new)
    plt.show()

    popt,_=curve_fit(gaussian_fit,logQ_new,M_new,p0=[7.0,0.1])

    return popt

def M_stats_cdf(logQ_array,M):
    
    M_cdf=_cummulative_distribution(logQ_array,mulitplicity=M)
    
    M_values=numpy.array([x[1] for x in M_cdf])
    logQ_mean=logQ_array[find_nearest(M_values,0.5)]
    
    logQ_sigma_right=logQ_array[find_nearest(M_values,0.6875)]
    logQ_sigma_left=logQ_array[find_nearest(M_values,0.3125)]
    return logQ_mean,logQ_sigma_left,logQ_sigma_right    
    # print(logQ_mean,logQ_sigma)
    # M_normal=norm.pdf(logQ_array,loc=logQ_mean,scale=logQ_sigma)/max(norm.pdf(logQ_array,loc=logQ_mean,scale=logQ_sigma))
    # plt.plot(logQ_array,M_normal,label='false')
    # plt.plot(logQ_array,M,color='green',label='true')
    # plt.legend()
    # plt.show()

if __name__=='__main__':

    s=['1', '8', '12', '13', '17', '20', '25', '28', '32', '36', '39', '43', '44', '47', '48', '50', '54', '56', '67', '70', '73', '76', '79', '80', '81', '83', '84', '85', '86', '88', '92', '93', '94', '95', '96', '106', '109', '120', '123', '126', '137']
    method=sys.argv[1]
    x=numpy.linspace(5,12,100000)

    with  open('all_pdf_data.pickle','rb') as f:
        D=pickle.load(f)

    M=D['M']
    M_interp=interpolate.InterpolatedUnivariateSpline(x,M,k=1)
    M_Z=M_interp.integral(5,12)
    M=M/M_Z
    M_mean,M_sigma_left,M_sigma_right=M_stats_cdf(x,M)
    sigma=(M_sigma_right-M_sigma_left)/2
    print(M_mean,sigma)

    P=p_value(M_mean,M_interp,M_Z)
    
    f=open(f'{method}.txt','w',1)
    f.write('System\tp-value(s)\n')
    for system in s:
        print(f'\nSystem: {system}')
        value=P(system,method)
        print(f'Obtained values: {value}')
        if method=='p_E':
            f.write(system+'\t'+repr(value)+'\n')
        else:f.write(system+'\t'+repr(value)+'\n')

    f.close()

