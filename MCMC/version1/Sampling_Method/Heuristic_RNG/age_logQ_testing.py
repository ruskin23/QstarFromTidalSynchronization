import numpy
import matplotlib.pyplot as plt
from utils import cummulative_distribution
from scipy.stats.stats import pearsonr
import itertools


M=500

def norm(x,mu,sigma):
    return (1/sigma*numpy.sqrt(2*3.14))*numpy.exp(-0.5*((x-mu)/(sigma))**2)

def sample(mu,
           sigma):

    return numpy.random.normal(loc=mu,scale=sigma)


def RMSE(y1,y2,desired_corr):

    corr=numpy.corrcoef(y1,y2)[0][1]

    return abs(corr-desired_corr)


def Algo(y1,y2,desired_corr,RMSE0):

    T=0.01
    dT=0.9
    gamma=0.01

    RMSE_i=RMSE0
    y2_best=numpy.copy(y2)
    y2_best_prev=numpy.copy(y2)

    N=0
    i=0

    while i<5:
        n_total=0
        n_accepted=0
        for k in range(5*M):
            N=N+1
            n_total=n_total+1

            i1=numpy.random.randint(0,M)
            i2=numpy.random.randint(0,M)

            y2[i1],y2[i2]=y2[i2],y2[i1]
            #import pdb; pdb.set_trace()
            RMSE_new=RMSE(y1,y2,desired_corr)

            delta=RMSE_i-RMSE_new

            with open('analysis.txt','a',1) as f:
                f.write(repr(i)+'\t'+
                        repr(N)+'\t'+
                        repr(RMSE_i)+'\t'+
                        repr(RMSE_new)+'\n')

            if delta>0 or delta==0:
                RMSE_i=RMSE_new
                n_accepted=n_accepted+1
            else:
                U=numpy.random.random(1)[0]
                P=numpy.exp(delta/T)
                if P>U:
                    RMSE_i=RMSE_new
                    n_accepted=n_accepted+1
                else:
                    y2[i2],y2[i1]=y2[i1],y2[i2]


        if RMSE(y1,y2,desired_corr)<RMSE(y1,y2_best,desired_corr):
            y2_best=numpy.copy(y2)

        if RMSE(y1,y2_best,desired_corr)<RMSE(y1,y2_best_prev,desired_corr):
            i=0

        r=n_accepted/n_total
        if r<gamma:
            i=i+1

        y2_best_prev=numpy.copy(y2_best)

        T=dT*T

    return RMSE_i

logQ_mean=9.037813422896482
logQ_sigma=0.5

logQ_samples=[]
logQ_prob=[]

for i in range(M):
    s=sample(logQ_mean,logQ_sigma)
    p=norm(s,logQ_mean,logQ_sigma)
    logQ_samples=numpy.append(logQ_samples,s)
    logQ_prob=numpy.append(logQ_prob,p)


age_mean=1.0491326397223144
age_sigma=1.0

age_samples=[]
age_prob=[]

age_samplefile='../../../mcmc_mass_age/samples/updated_samples/MassAgeFehSamples_137.txt'

N=0
MM=[]
A=[]

mass_mean=1.0868693079612135
mass_sigma=1.0

feh_mean=-0.041740014503029985
feh_sigma=0.3

with open(age_samplefile,'r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        mass=float(x[0])
        age=float(x[1])
        feh=float(x[2])
        multiplicity=float(x[3])
        modified_multiplicity=multiplicity*numpy.exp(-0.5*((age-age_mean)/(age_sigma))**2)*numpy.exp(-0.5*((mass-mass_mean)/(mass_sigma))**2)*numpy.exp(-0.5*((feh-feh_mean)/(feh_sigma))**2)
        N=N+modified_multiplicity

        A=numpy.append(A,age)
        MM=numpy.append(MM,modified_multiplicity)



MM=MM/N
for i in range(M):
    U=numpy.random.random(1)[0]
    for k in range(len(A)):
        if MM[k]>U:
            age_samples=numpy.append(age_samples,A[k])
            age_prob=numpy.append(age_prob,MM[k])
            break
        else:U=U-MM[k]



age_old_samples=numpy.zeros(len(age_samples))
for i in range(len(age_samples)):
    age_old_samples[i]=age_samples[i]
print(numpy.corrcoef(logQ_samples,age_samples)[0][1])

desired_corr=0.13961937
RMSE0=RMSE(logQ_samples,age_samples,desired_corr)

RMSE_f=Algo(logQ_samples,age_samples,desired_corr,RMSE0)

print(numpy.corrcoef(logQ_samples,age_samples)[0][1])
print(RMSE_f)

cdf_old=cummulative_distribution(age_old_samples,age_prob)()
cdf_new=cummulative_distribution(age_samples,age_prob)()

for co,cn in zip(cdf_old,cdf_new):
    print('{}\t{}'.format(co,cn))

plt.plot(*zip(*cdf_old),color='r')
plt.plot(*zip(*cdf_new),color='g')
plt.show()
