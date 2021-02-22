import numpy 
import pickle
import matplotlib.pyplot as plt
from scipy import integrate
from scipy import interpolate
from utils import _get_filename,_get_chain,_fill_parameters,_cummulative_distribution

def cummulative_distribution(x,y):

    P=0
    cdf=[]
    for i in range(100000):
        P=P+y[i]
        cdf.append(P)
    cdf=cdf/max(cdf)

    for i,c in enumerate(cdf):
        if abs(c-0.5)<1e-4:
            return x[i]
    

def get_chain(system):
    CHAIN=numpy.zeros(1)

    for c in ['ganymede','stampede']:
        for i in range(1,6):
            filename=_get_filename(system,c,i)
            filled_file=_fill_parameters(filename)
            chain=_get_chain('logQ',filled_file)
            if filename!=f'ganymede/MCMC_{system}/accepted_parameters_1.txt':
                start_point=chain[0]
                count_repeat=numpy.count_nonzero(chain==start_point)
                chain=chain[count_repeat:]
            CHAIN=numpy.concatenate((CHAIN,chain),axis=None)

    return CHAIN


def kernel_expan(x,x_i,h):
    y=(x-x_i)/h
    y=0.75*(1-y**2)
    y[y<0]=0
    return y

def kernel_gauss(x,x_i,h):

    y=(x-x_i)/h
    return (1.0/(numpy.sqrt(2*numpy.pi)))*numpy.exp(-(y**2)/2) 

def kde(x,x_array,h=None):

    n=len(x_array)
    if h is None:h=3.5*numpy.power(n,-1.0/3)

    f=0
    for i in range(n):
        f=f+kernel_gauss(x,x_array[i],h)
    return f/(n*h)
    

def p_value_E(c,logQ_samples,logQ_cdf,logQ_median):
    
    # lower_limit=kde(5.0,logQ_samples)
    # upper_limit=kde(12.0,logQ_samples)

    # if numpy.logical_or(numpy.logical_and(lower_limit==0,
    #                                     upper_limit==0),
    #                 numpy.logical_and(lower_limit!=0,
    #                                     upper_limit!=0)
    #                 ):

    x=numpy.linspace(5,12,100000)
    c_quantile=full_cdf_logQ(c, logQ_samples, logQ_cdf)
    other_quantile=1-c_quantile

    for v in x:
        q=full_cdf_logQ(v,logQ_samples,logQ_cdf)
        if abs(q-other_quantile)<1e-3:break

    p=[]
    if c>logQ_median:
        p.append(integrate.quad(kde,5,v,args=(logQ_samples))[0])
        p.append(integrate.quad(kde,c,12,args=(logQ_samples))[0])
    else:
        p.append(integrate.quad(kde,5,c,args=(logQ_samples))[0])
        p.append(integrate.quad(kde,v,12,args=(logQ_samples))[0])

    
    return p

def full_cdf_logQ(x,logQ_samples,cdf_logQ):

    logQ_interp=interpolate.interp1d(logQ_samples,cdf_logQ)

    if x<min(logQ_samples):return 0
    elif x>max(logQ_samples):return 1
    else:return logQ_interp(x)    


def integrand_E_p(x,logQ_samples,cdf_logQ,interp_M):

    cdf_p=full_cdf_logQ(x,logQ_samples,cdf_logQ)
    q=min(cdf_p,1-cdf_p)
    print(f'x={x};q={q*interp_M(x)}')
    return q*interp_M(x)

def E_p_value(logQ_samples,cdf_logQ,pdf_M):

    logQ_array=numpy.linspace(5,12,100000)
    interp_M=interpolate.interp1d(logQ_array,pdf_M)
    p=integrate.quad(integrand_E_p,5,12,args=(logQ_samples,cdf_logQ,interp_M))

    return p[0]


if __name__=='__main__':
    with  open('all_pdf_data.pickle','rb') as f:
        D=pickle.load(f)

    s=['85', '73', '76', '96', '92', '81', '80', '36', '93', '83', '84', '94', '32', '79', '106', '123', '50', '47', '39', '56', '126', '54', '109', '44', '48', '17', '70', '8', '12', '88', '67', '20', '95', '25', '137', '120', '86', '43', '28', '13']
    M=D['M']/max(D['M'])
    logQ_array=numpy.linspace(5,12,100000)

    M_mean=cummulative_distribution(logQ_array,M)
    # f=open('E_p_table.txt','w',1)
    f=open('p_E_table_temp.txt','w',1)
    f.write('System\tp-value\n')
    # p_E=[]
    # E_p=[]
    # system_names=[]
    for system in s:
    #     with open('p_E_table.txt','r') as f:
    #         for lines in f:
    #             x=lines.split()
    #             if x[0]==system:
    #                 system_names.append(int(x[0]))
    #                 p_E.append(float(x[1]))
    #                 break

    #     with open('E_p_table.txt','r') as f:
    #         for lines in f:
    #             x=lines.split()
    #             if x[0]==system:
    #                 E_p.append(float(x[1]))
    #                 break

        print(f'System = {system}')
        logQ_raw=get_chain(system)
        cdf=_cummulative_distribution(logQ_raw)
        Z=list(zip(*cdf))
        logQ_values=list(Z[0])
        P=list(Z[1])
        logQ=[]
        P_idx=[]
        for q in logQ_values:
            if q not in logQ:logQ.append(q)
            else:
                print(q)
                P_idx.append(logQ.index(q))
        logQ_samples=numpy.array(logQ)
        logQ_median=numpy.median(logQ_samples)
        P=numpy.array(P)
        P=numpy.delete(P,P_idx)

        p=p_value_E(M_mean,logQ_samples,P,logQ_median)
        # p=E_p_value(logQ_samples,P,M)
        f.write(system+'\t'+repr(p[0])+'\t'+repr(p[1])+'\n')

    # plt.scatter(system_names,p_E,color='red',label='p_E')
    # plt.scatter(system_names,E_p,color='green',label='E_p')
    # plt.legend()
    # plt.show()
    # f.close()