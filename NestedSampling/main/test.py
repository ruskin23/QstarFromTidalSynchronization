import numpy
import scipy
from scipy import stats
import matplotlib.pyplot as plt

def ptform(u):

    x=numpy.array(u)

    for i,s in enumerate(sampling_parameters):

        if s[-1]=='Normal':
            mean=s[1]
            sigma=s[2]
            x[i]=scipy.stats.norm.ppf(u[i],loc=mean,scale=sigma)
        elif s[-1]=='Turncated_Normal':
            mean=s[1]
            sigma=s[2]
            low=(s[3]-mean)/sigma
            high=(s[4]-mean)/sigma
            x[i]=scipy.stats.truncnorm.ppf(u[i],low,high,loc=mean,scale=sigma)
        elif s[-1]=='Uniform':
            x[i]=(s[2]-s[1])*u[i] + s[1]

    return x

catalog_file='SpinlogQCatalog_el0.4.txt'

with open(catalog_file,'r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        if x[0]=='39':
            teff_value=float(x[2])
            teff_error=float(x[3])
            feh_value=float(x[4])
            feh_error=float(x[5])
            logg_value=float(x[10])
            logg_error=float(x[15])
            Porb_value=float(x[6])
            Porb_error=float(x[7])
            eccentricity_value=float(x[8])
            eccentricity_error=float(x[9])
            Pspin_value=float(x[12])
            Pspin_error=float(x[13])
            mass_ratio=float(x[14])
            break


sampling_parameters = [('Porb',Porb_value,Porb_error,'Normal'),
                        ('feh',feh_value,feh_error,-1.014,0.537,'Turncated_Normal'),
                        ('eccentricity',eccentricity_value,eccentricity_error,0.0,0.45,'Turncated_Normal'),
                        ('Wdisk',2*scipy.pi/14,2*scipy.pi/1.4,'Uniform'),
                        ('logQ',5.0,12.0,'Uniform'),
                        ('primary_mass',0.5,1.2,'Uniform'),
                        ('age',1e-3,10.0,'Uniform')]


# r=numpy.random.rand(10,7)
# p=numpy.array(list(map(ptform,r)))
# print(p)
# print('break')
# outfile='temp'
# for i in range(len(p)):
#     numpy.savez(outfile,p=p,pi=p[i])

# npzfile=numpy.load('temp.npz',allow_pickle=True)
# print(npzfile['pi'])

#sub-linear scaling
K=numpy.linspace(50,1500,1000)
M=numpy.array([10,20,30,40,50])
color=['r','g','b','y','k']

for m,c in zip(M,color):
    S=K*numpy.log(1+(m/K))
    plt.plot(S/K,m/K,color=c,label=str(m))
plt.legend()
plt.show()