import numpy
import matplotlib.pyplot as plt
from utils import cummulative_distribution
from scipy.stats.stats import pearsonr



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

for i in range(M):
    y1=numpy.append(y1,sample(mu1,sigma1))
    y2=numpy.append(y2,sample(mu2,sigma2))

print('Old Corr = ',numpy.corrcoef(y1,y2)[0][1])
RMSE_i=RMSE(y1,y2)
print('Old RMSE = ',RMSE_i)

T=0.01
dT=0.9
gamma=0.01

y2_best=y2
y2_best_prev=y2

i=0
while i<5:
    n_total=0
    n_accepted=0
    for k in range(5*M):

        n_total=n_total+1

        i1=numpy.random.randint(0,100)
        i2=numpy.random.randint(0,100)

        y2[i1],y2[i2]=y2[i2],y2[i1]

        RMSE_new=RMSE(y1,y2)

        delta=RMSE_i-RMSE_new

        if delta>0:
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



