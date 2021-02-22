import sys
import numpy 
import pickle
import matplotlib.pyplot as plt
from scipy import integrate
from scipy import interpolate
from scipy.stats import norm
from utils import _get_filename,_get_chain,_fill_parameters,_cummulative_distribution



class p_value():

    def __init__(self,M_mean,M_interp):

        
        self.clusters=['ganymede','stampede']

        self.logQ_array=numpy.linspace(5,12,100000)
        self.two_sided=False
        self.M_mean=M_mean
        self.M_interp=M_interp


    def get_chain(self,system):

        self.CHAIN=numpy.zeros(1)
        
        for c in self.clusters:
            for i in range(1,6):
                filename=_get_filename(system,c,i)
                filled_file=_fill_parameters(filename)
                chain=_get_chain('logQ',filled_file)
                if filename!=f'ganymede/MCMC_{system}/accepted_parameters_1.txt':
                    start_point=chain[0]
                    count_repeat=numpy.count_nonzero(chain==start_point)
                    chain=chain[count_repeat:]
                self.CHAIN=numpy.concatenate((self.CHAIN,chain),axis=None)
        self.CHAIN=self.CHAIN[1:]
    
    def get_samples(self,system):

        self.get_chain(system)
        cdf=_cummulative_distribution(self.CHAIN)
        Z=list(zip(*cdf))
        self.logQ_values=numpy.array(list(Z[0]))
        self.logQ_cdf=numpy.array(list(Z[1]))
        self.logQ_median=numpy.median(self.logQ_values)
        self.logQ_interp=interpolate.interp1d(self.logQ_values,self.logQ_cdf)
        

    def get_cdf_statistic(self,quantile=None,value=None):
        if quantile is not None:
            for i,c in enumerate(self.logQ_cdf):
                if abs(c-quantile)<1e-3:return self.logQ_values[i]
        
        if value is not None:
            if value<min(self.logQ_values):return 0
            elif value>max(self.logQ_values): return 1
            else:return self.logQ_interp(value)
    
    def logQ_pdf(self,x,h=None):

        N=len(self.logQ_values)
        if h is None:h=3.5*numpy.power(N,-1.0/3)

        f_h=0
        for i in range(N):
            y=(x-self.logQ_values[i])/h
            f_h=f_h+norm.pdf(y)

        return f_h/(N*h)

    def p_E(self):

        M_quantile=self.get_cdf_statistic(value=self.M_mean)
        print(f'M_mean={self.M_mean} quantile={M_quantile}')
        M_other_quantile=1-M_quantile
        q_value=self.get_cdf_statistic(quantile=M_other_quantile)
        print(f'q-value={q_value} quantile={M_other_quantile}')

        lower_limit=self.logQ_pdf(5.0)
        upper_limit=self.logQ_pdf(12.0)

        print(f'Lower Limit = {lower_limit}; Upper Limit = {upper_limit}')

        if numpy.logical_or(numpy.logical_and(lower_limit<1e-5,
                                            upper_limit<1e-5),
                        numpy.logical_and(lower_limit>1e-3,
                                            upper_limit>1e-3)
                        ):self.two_sided=True
        else:self.two_sided==False

        print(f'Two-sided={self.two_sided}')
        p=[]
        if self.M_mean>self.logQ_median:
            print('Integrating for first condition')
            if self.two_sided==True: p.append(integrate.quad(self.logQ_pdf,5,q_value)[0])
            p.append(integrate.quad(self.logQ_pdf,self.M_mean,12)[0])
        else:
            print('Integrating for second condition')
            p.append(integrate.quad(self.logQ_pdf,5,self.M_mean)[0])
            if self.two_sided==True: p.append(integrate.quad(self.logQ_pdf,q_value,12)[0])
        print('Done')
        return p

    def integrand_E_p(self,x):

        cdf_p=self.get_cdf_statistic(value=x)
        q=min(cdf_p,1-cdf_p)
        return q*self.M_interp(x)

    def E_p(self):

        p=integrate.quad(self.integrand_E_p,5,12)
        return p[0]


    def __call__(self,system,method):

        
        self.get_samples(system)
        if method=='p_E':return self.p_E()
        elif method=='E_p':return self.E_p()
        else: raise ValueError


if __name__=='__main__':

    s=['85', '73', '76', '96', '92', '81', '80', '36', '93', '83', '84', '94', '32', '79', '106', '123', '50', '47', '39', '56', '126', '54', '109', '44', '48', '17', '70', '8', '12', '88', '67', '20', '95', '25', '137', '120', '86', '43', '28', '13']
    method=sys.argv[1]

    with  open('all_pdf_data.pickle','rb') as f:
        D=pickle.load(f)

    M=D['M']/max(D['M'])
    logQ_array=numpy.linspace(5,12,100000)
    for i,m in enumerate(M):
        if abs(m-1.0)<1e-4:M_mean=logQ_array[i]

    M_interp=interpolate.interp1d(logQ_array,M)

    P=p_value(M_mean,M_interp)

    f=open(f'{method}.txt','w',1)
    f.write('System\tp-value(s)\n')
    for system in s:
        print(f'\nSystem: {system}')
        value=P(system,method)
        print(f'Obtained values: {value}')
        if method=='p_E':
            if len(value)==1:f.write(system+'\t'+repr(value[0])+'\n')
            else:f.write(system+'\t'+repr(value[0])+'\t'+repr(value[1])+'\t'+repr(P.two_sided)+'\n')
        else:f.write(system+'\t'+repr(value)+'\n')

    f.close()

