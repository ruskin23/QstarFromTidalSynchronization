#Testing for bivariate normal

import numpy
import matplotlib.pyplot as plt
from utils import cummulative_distribution
from scipy.stats.stats import pearsonr

def norm(x,mu,sigma):
    return (1/sigma*numpy.sqrt(2*3.14))*numpy.exp(-0.5*((x-mu)/(sigma))**2)

def sample(mu,
           sigma):

    return numpy.random.normal(loc=mu,scale=sigma)


def RMSE(y1,y2):

    corr=numpy.corrcoef(y1,y2)[0][1]
    desired_corr=0.8

    return abs(corr-desired_corr)

mu1=3.0
mu2=7.7

sigma1=0.04
sigma2=0.08

M=100

y1=[]
y2=[]
y2_prob=[]


for i in range(M):
    y1=numpy.append(y1,sample(mu1,sigma1))
    y2_sample=sample(mu2,sigma2)
    y2=numpy.append(y2,y2_sample)
    y2_prob=numpy.append(y2_prob,y2_sample)

cdf=cummulative_distribution(y2,y2_prob)()
plt.plot(*zip(*cdf))

index=numpy.random.randint(0,M)
print(y1[index])
print(y2[index])

print('Old Corr = ',numpy.corrcoef(y1,y2)[0][1])
RMSE_i=RMSE(y1,y2)
print('Old RMSE = ',RMSE_i)

T=0.01
dT=0.9
gamma=0.01

y2_best=y2
y2_best_prev=y2

i=0


with open('analysis.txt','w',1) as f:
    f.write('i'+'\t'+'N'+'\t'+'delta'+'\n')

N=0
while i<5:
    n_total=0
    n_accepted=0
    for k in range(5*M):
        N=N+1
        n_total=n_total+1

        i1=numpy.random.randint(0,M)
        i2=numpy.random.randint(0,M)

        y2[i1],y2[i2]=y2[i2],y2[i1]

        RMSE_new=RMSE(y1,y2)

        delta=RMSE_i-RMSE_new

        with open('analysis.txt','a',1) as f:
            f.write(repr(i)+'\t'+
                    repr(N)+'\t'+
                    repr(numpy.log10(RMSE_i))+'\n')

        if delta>0 or delta==0:
            RMSE_i=RMSE_new
            n_accepted=n_accepted+1
            continue
        else:
            U=numpy.random.random(1)[0]
            P=numpy.exp(delta/T)
            if P>U:
                RMSE_i=RMSE_new
                n_accepted=n_accepted+1
                continue
            else:
                y2[i2],y2[i1]=y2[i1],y2[i2]


    if RMSE(y1,y2)<RMSE(y1,y2_best):
        y2_best=y2

    if RMSE(y1,y2_best)<RMSE(y1,y2_best_prev):
        i=0

    r=n_accepted/n_total
    if r<gamma:
        i=i+1

    y2_best_prev=y2_best

    T=dT*T

print('New Corr = ',numpy.corrcoef(y1,y2)[0][1])
print('New RMSE = ',RMSE_i)

print(y1[index])
print(y2[index])

print(numpy.cov(y1,y2))

cdf=cummulative_distribution(y2,y2_prob)()
plt.plot(*zip(*cdf))
plt.show()
