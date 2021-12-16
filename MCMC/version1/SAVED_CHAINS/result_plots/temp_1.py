import  sys
import numpy
import pickle
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.stats import norm
from utils import _get_filename,_get_chain,_fill_parameters,_cummulative_distribution



# x=numpy.linspace(5,12,100000)
# M_interp=interpolate.interp1d(x,M)

# plt.plot(x,M)
# plt.plot(x,M_interp(x))
# plt.show()
def get_chain(system):

    CHAIN=[]
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
    CHAIN=CHAIN[1:]
    return CHAIN

def get_samples(system):

    CHAIN=get_chain(system)
    cdf=_cummulative_distribution(CHAIN)
    Z=list(zip(*cdf))
    logQ_values=numpy.array(list(Z[0]))
    logQ_cdf=numpy.array(list(Z[1]))
    return logQ_values,logQ_cdf


def get_quantile(v,logQ_values,logQ_cdf_interp):
    if v<min(logQ_values):return 0
    elif v>max(logQ_values):return 1
    else: return logQ_cdf_interp(v)
    
def get_pdf(x,h,N):

    f_h=0
    for i in range(N):
        y=(x-logQ_values[i])/h
        f_h=f_h+norm.pdf(y)

    return f_h/(N*h)

def find_nearest(array, value):
    idx = (numpy.abs(array - value)).argmin()
    return idx,array[idx]


with open('all_pdf_data.pickle','rb') as f:
    D=pickle.load(f)
M=D['M']
M=M/max(M)
X=numpy.linspace(5,12,100000)
for i,m in enumerate(M):
    if abs(m-1.0)<1e-4:M_mean=X[i]

x=numpy.linspace(5,12,100000)

s=['85', '73', '76', '96', '92', '81', '80', '36', '93', '83', '84', '94', '32', '79', '106', '123', '50', '47', '39', '56', '126', '54', '109', '44', '48', '17', '70', '8', '12', '88', '67', '20', '95', '25', '137', '120', '86', '43', '28', '13']

s=numpy.reshape(numpy.array(s),(8,5))
fig,axs=plt.subplots(8,5,sharex='col',sharey='row',gridspec_kw={'hspace': 0, 'wspace': 0})


for i in range(8):
    for j in range(5):
        system=s[i][j]
        print(system)

        logQ_values,logQ_cdf=get_samples(system)

        N=len(logQ_values)
        h=3.5*numpy.power(N,-1.0/3)
        f_h = get_pdf(x,h,N)/(N*h)
       

        with open('p_E.txt','r') as f:
            for lines in f:
                l=lines.split()
                if l[0]==system:
                    p_left=float(l[1])
                    p_right=float(l[2])          
                    x_lower=float(l[4])
                    x_upper=float(l[6])
                    break

        x_lower_idx,x_lower_value=find_nearest(x,x_lower)
        x_upper_idx,x_upper_value=find_nearest(x,x_upper)

        print(x_lower_value,x_upper_value)

        x_lower_array=x[0:x_lower_idx]
        x_upper_array=x[x_upper_idx:-1]

        y_lower_array=f_h[0:x_lower_idx]
        y_upper_array=f_h[x_upper_idx:-1]

        if M_mean>numpy.median(logQ_values):
            color_left='blue'
            color_right='green'
        else:
            color_left='green'
            color_right='blue'

        axs[i,j].plot(x,f_h,color='black',label=f'Sytem_{system}_pdf')
        axs[i,j].plot(X,M,linestyle='--',color='darkblue',label='combined_pdf')
        axs[i,j].fill_between(x_lower_array,y_lower_array,color=color_left,label=f'p={round(p_left,5)}')
        axs[i,j].fill_between(x_upper_array,y_upper_array,color=color_right,label=f'p={round(p_right,5)}')
        axs[i,j].grid(True)
        # axs[i,j].legend()
        axs[i,j].tick_params(axis='x', labelsize= 3)
        axs[i,j].tick_params(axis='y', labelsize= 3)


for ax in axs.flat:
    ax.label_outer()


plt.savefig('pdf/p_values/all_pdf_mean.png')

# for system in s:
#     print(system)

#     logQ_values,logQ_cdf=get_samples(system)
#     logQ_cdf_interp=interpolate.interp1d(logQ_values,logQ_cdf)

#     x=numpy.linspace(5,12,10000)
#     y=numpy.zeros(len(x))
#     for i in range(len(x)):
#         y[i]=get_quantile(x[i],logQ_values,logQ_cdf_interp)

#     logQ_values_interp=interpolate.interp1d(logQ_cdf,logQ_values)
#     y1=numpy.linspace(min(logQ_cdf),max(logQ_cdf),1000)
#     x1=logQ_values_interp(y1)


#     plt.plot(x,y)
#     plt.scatter(x1,y1)
#     plt.show()

#     x=numpy.linspace(5,12,50000)
#     N=len(logQ_values)
#     h=3.5*numpy.power(N,-1.0/3)
#     f_h = get_pdf(x,h,N)/(N*h)
#     f_h=f_h/max(f_h)

#     with open('p_E.txt','r') as f:
#         for lines in f:
#             l=lines.split()
#             if l[0]==system:
#                 p_left=float(l[1])
#                 p_right=float(l[2])          
#                 x_lower=float(l[4])
#                 x_upper=float(l[6])

#     x_lower_idx,x_lower_value=find_nearest(x,x_lower)
#     x_upper_idx,x_upper_value=find_nearest(x,x_upper)

#     print(x_lower_value,x_upper_value)

#     x_lower_array=x[0:x_lower_idx]
#     x_upper_array=x[x_upper_idx:-1]

#     y_lower_array=f_h[0:x_lower_idx]
#     y_upper_array=f_h[x_upper_idx:-1]

#     if M_mean>numpy.median(logQ_values):
#         color_left='blue'
#         color_right='green'
#     else:
#         color_left='green'
#         color_right='blue'

#     plt.plot(x,f_h,color='black',label=f'Sytem_{system}_pdf')
#     plt.plot(X,M,linestyle='--',color='darkblue',label='combined_pdf')
#     plt.fill_between(x_lower_array,y_lower_array,color=color_left,label=f'p={round(p_left,5)}')
#     plt.fill_between(x_upper_array,y_upper_array,color=color_right,label=f'p={round(p_right,5)}')
#     plt.grid(True)
#     plt.legend(loc='upper right')
#     plt.savefig(f'pdf/p_values/System_{system}.png')
#     plt.close()