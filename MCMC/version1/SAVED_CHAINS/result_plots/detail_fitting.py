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
    


def get_pdfs(s,save_plots=False,dump_pdf=False):

    M=numpy.ones(100000)
    PDF=[]
    for system in s:
        print(system)
        logQ_raw=get_chain(system)

        x=numpy.linspace(5,12,100000)
        f=kde(x,logQ_raw)
        PDF.append(f)
        M=M*f
        if save_plots==True:
            plt.xlim(5,12)
            plt.plot(x,f,color='blue')
            plt.savefig(f'pdf/gauss_kde_new/system_{system}.png')
            plt.close()
    
    dataset=dict()
    for name,array in zip(s,PDF):
        dataset[name]=array
    dataset['M']=M

    if dump_pdf==True:
        with open('all_pdf_data.pickle','wb') as f:
            pickle.dump(dataset,f)

    return dataset        



def discard_low_e_system(s):
    
    Q=numpy.linspace(5,12,100000)

    s_remove=[]
    with open('../SpinlogQCatalog_el0.4.txt','r') as f:
        next(f)
        for lines in f:
            x=lines.split()
            if float(x[8])<0.1 and x[0] not in ['106','86']:
                s_remove.append(x[0])
    
    s_new=[systems for systems in s if systems not in s_remove]

    with open('all_pdf_data.pickle','rb') as f:
        D=pickle.load(f)
    
    M_old=D['M']/interpolate.InterpolatedUnivariateSpline(Q,D['M']).integral(5,12)
    M_new=numpy.ones(100000)
    for system in s_new:
        M_new=M_new*D[system]
    M_new=M_new/interpolate.InterpolatedUnivariateSpline(Q,M_new).integral(5,12)

    plt.plot(Q,M_old,color='g')
    plt.plot(Q,M_new,color='r')
    plt.show()
    

def plot_parameters(system,param1,param2):

    with open('../complete_chains.pickle','rb') as f:
        D=pickle.load(f)
    
    for s,_ in D.items():
        if s==system:
            p_1=D[s][param1]
            p_2=D[s][param2]
            break


    plt.scatter(p_1,p_2)
    # p_1_min,p_2_min=min(p_1),min(p_1)
    # p_1_max,p_2_max=max(p_1),max(p_2)

    # p_1_bins=numpy.linspace(p_1_min,p_1_max,100)
    # p_2_bins=numpy.linspace(p_2_min,p_2_max,100)

    # plt.hist2d(p_1,p_2,bins=[p_1_bins,p_2_bins])
    # plt.ylim(8.0,15.0)

    # heatmap, xedges, yedges = numpy.histogram2d(p_1, p_2, bins=50)
    # print(xedges,yedges)
    # extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

    # plt.clf()
    # plt.imshow(heatmap.T, extent=extent, origin='lower')


    # spin_min=10.778-0.38400000000000034
    # spin_max=10.778+0.38400000000000034
    # spin_value=10.778  

    # plt.axhline(y=spin_min,color='r',linestyle='--')
    # plt.axhline(y=spin_max,color='r',linestyle='--')
    # plt.axhline(y=spin_value,color='k')
    plt.show()

if __name__=='__main__':

    s=['1', '8', '12', '13', '17', '20', '25', '28', '32', '36', '39', '43', '44', '47', '48', '50', '54', '56', '67', '70', '73', '76', '79', '80', '81', '83', '84', '85', '86', '88', '92', '93', '94', '95', '96', '106', '109', '120', '123', '126', '137']


    # discard_low_e_system(s)
    plot_parameters('1','logQ','Spin')
