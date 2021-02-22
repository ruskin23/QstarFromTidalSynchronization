import scipy
from scipy import interpolate
from scipy.misc import derivative
from utils import _get_filename,_get_chain,_fill_parameters,_cummulative_distribution
import numpy
import matplotlib.pyplot as plt
import sys
from collections import OrderedDict 


def get_chain(system):
    CHAIN=numpy.zeros(1)

    for c in ['ganymede','stampede']:
        for i in range(1,6):
            filename=_get_filename(system,c,i)
            filled_file=_fill_parameters(filename)
            chain=_get_chain('logQ',filled_file)
            if i!=1:
                start_point=chain[0]
                count_repeat=numpy.count_nonzero(chain==start_point)
                chain=chain[count_repeat:]
            CHAIN=numpy.concatenate((CHAIN,chain),axis=None)

    return CHAIN[1:]




# s=['85', '73', '76', '96', '92', '81', '80', '36', '93', '83', '84', '94', '32', '79', '106', '123', '50', '47', '39', '56', '126', '54', '109', '44', '48', '17', '70', '8', '12', '88', '67', '20', '95', '25', '137', '120', '86', '43', '28', '13']
# for system in s:
system=sys.argv[1]
logQ_values=get_chain(system)
cdf=_cummulative_distribution(logQ_values)
Z=list(zip(*cdf))
logQ=list(Z[0])
P=list(Z[1])
logQ_new=[]
for q in logQ:
    if q not in logQ_new:logQ_new.append(q)
    else:
        print(q)
        del P[logQ.index(q)]
if (all(i < j for i, j in zip(logQ_new, logQ_new[1:]))):
    print("Yes, List is sorted.")
else:
    print("No, List is not sorted.")
print(f'system={system} logQ_min={min(logQ_new)} logQ_max={max(logQ_new)}')
plt.plot(logQ_new,P,color='blue')

# box=numpy.ones(5)/5
# P_smooth=numpy.convolve(P,box,mode='same')
# plt.plot(logQ_new,P_smooth)


spl=interpolate.interp1d(logQ_new,P)
x=numpy.linspace(min(logQ_new),max(logQ_new),10000)
y=spl(x)
plt.plot(x,y,color='green')
# plt.xlim(5.4,5.5)
# plt.show()


deriv=[]
x_new=[]
for i in range(1,len(x)):
    dx=x[i]-x[i-1]
    x_new.append(x[i-1]+dx/2)
    deriv.append( (y[i]-y[i-1]) / dx )
deriv=deriv/max(deriv)
plt.plot(x_new,deriv)
# plt.xlim(5.4,5.5)
plt.ylim(0,1)
plt.show()

    # max_idx=numpy.where(deriv==numpy.amax(deriv))
    # deriv_gauss=deriv[0:2*max_idx[0][0]]
    # x_gauss=x_new[0:2*max_idx[0][0]]
    # plt.scatter(x_gauss,deriv_gauss)
    # plt.show()

    # deriv_tail=deriv[2*max_idx[0][0]:-1]
    # x_tail=x_new[2*max_idx[0][0]:-1]
    # plt.scatter(x_tail,deriv_tail)
    # plt.show()

# n,bins,patches=plt.hist(logQ,bins='auto',alpha=0.7, rwidth=0.85)
    # plt.ylim(200,1600)
# plt.show()