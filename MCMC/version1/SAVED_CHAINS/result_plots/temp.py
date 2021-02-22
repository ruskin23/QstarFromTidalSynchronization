import numpy
import pickle
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.stats import norm
from utils import _get_filename,_get_chain,_fill_parameters,_cummulative_distribution


with open('all_pdf_data.pickle','rb') as f:
    D=pickle.load(f)
# M=D['M']
# M=M/max(M)
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
    return logQ_values

logQ_values=get_samples('73')

# x=numpy.linspace(5,12,1000)

x=numpy.linspace(8,12,10)
N=len(logQ_values)
h=3.5*numpy.power(N,-1.0/3)
f_h=0
for i in range(N):
    y=(x-logQ_values[i])/h
    f_h=f_h+norm.pdf(y)

f_h = f_h/(N*h)
print(f_h)
# plt.plot(x,f_h)
# plt.show()