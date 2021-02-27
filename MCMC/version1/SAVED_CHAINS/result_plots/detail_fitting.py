import scipy
from scipy import interpolate
from scipy.misc import derivative
from utils import _get_filename,_get_chain,_fill_parameters,_cummulative_distribution,adjust_chain
import numpy
import matplotlib.pyplot as plt
import sys
from collections import OrderedDict 
import pickle

def get_chain(system):
    # CHAIN=numpy.zeros(1)

    # for c in ['ganymede','stampede']:
    #     for i in range(1,6):
    #         filename=_get_filename(system,c,i)
    #         filled_file=_fill_parameters(filename)
    #         chain=_get_chain('logQ',filled_file)
    #         if i!=1:
    #             start_point=chain[0]
    #             count_repeat=numpy.count_nonzero(chain==start_point)
    #             chain=chain[count_repeat:]
    #         CHAIN=numpy.concatenate((CHAIN,chain),axis=None)
    # CHAIN=adjust_chain(system,CHAIN[1:])
    # return CHAIN


    with open('../complete_chains.pickle','rb') as f:
        D=pickle.load(f)
    for system_name,params in D.items():
        if system_name==system:
            for name,values in params.items():
                if name=='logQ':
                    return values

def kernel_expan(x,x_i,h):
    y=(x-x_i)/h
    y=0.75*(1-y**2)
    y[y<0]=0
    return y

def kernel_gauss(x,x_i,h):

    y=(x-x_i)/h
    return (1/(numpy.sqrt(2*numpy.pi)))*numpy.exp(-(y**2)/2) 

def kde(x,x_array,h=None):

    n=len(x_array)
    if h is None:h=3.5*numpy.power(n,-1.0/3)

    f=0
    for i in range(n):
        f=f+kernel_gauss(x,x_array[i],h)
    return f/(n*h)
    

s=['1', '8', '12', '13', '17', '20', '25', '28', '32', '36', '39', '43', '44', '47', '48', '50', '54', '56', '67', '70', '73', '76', '79', '80', '81', '83', '84', '85', '86', '88', '92', '93', '94', '95', '96', '106', '109', '120', '123', '126', '137']
M=numpy.ones(100000)
PDF=[]
for system in s:
    print(system)
    logQ_raw=get_chain(system)

    x=numpy.linspace(5,12,100000)
    f=kde(x,logQ_raw)
    PDF.append(f)
    M=M*f
    plt.xlim(5,12)
    plt.plot(x,f,color='blue')
    plt.savefig(f'pdf/gauss_kde_new/system_{system}.png')
    plt.close()
        

dataset=dict()
for name,array in zip(s,PDF):
    dataset[name]=array
dataset['M']=M

with open('all_pdf_data.pickle','wb') as f:
    pickle.dump(dataset,f)



    # cdf=_cummulative_distribution(logQ_raw)
    # Z=list(zip(*cdf))
    # logQ_values=list(Z[0])
    # P=list(Z[1])
    # logQ=[]
    # P_idx=[]
    # for q in logQ_values:
    #     if q not in logQ:logQ.append(q)
    #     else:
    #         print(q)
    #         P_idx.append(logQ.index(q))
    # logQ=numpy.array(logQ)
    # P=numpy.array(P)
    # P=numpy.delete(P,P_idx)
    # print(f'system={system} logQ_min={min(logQ)} logQ_max={max(logQ)}')
    # plt.plot(logQ,P,color='red')

    # q_diff=[]
    # p_diff=[]
    # r=[]
    # for i in range(1,len(logQ)):
    #     q_diff.append(logQ[i]-logQ[i-1])
    #     p_diff.append(P[i]-P[i-1])
    #     r.append((P[i]-P[i-1])/(logQ[i]-logQ[i-1]))
    # N=1
    # count=len([i for i in r if i>N])
    # print(f'total={len(logQ)},greater than {N}={count}, remaining={len(logQ)-count}')

    # g_idx=[idx+1 for idx,val in enumerate(r) if val>N]+[idx for idx,val in enumerate(r) if val>N]
    # logQ_new=[val for idx,val in enumerate(logQ) if idx not in g_idx]
    # P_new=[val for idx,val in enumerate(P) if idx not in g_idx]
    # print(len(logQ),len(logQ_new))
    # plt.plot(logQ_new,P_new)
    
    # spl=interpolate.UnivariateSpline(logQ,P,s=0.02)
    # x=numpy.linspace(min(logQ),max(logQ),10000)
    # y=spl(x)
    # plt.plot(x,y,color='green')
    # plt.show()

    # deriv=[]
    # x_new=[]
    # for i in range(1,len(x)):
    #     dx=x[i]-x[i-1]
    #     x_new.append(x[i-1]+dx/2)
    #     deriv.append( (y[i]-y[i-1]) / dx )
    # deriv=deriv/max(deriv)
    # plt.plot(x_new,deriv,'black')
    # plt.show()
