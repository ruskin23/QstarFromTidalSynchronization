import sys
import numpy
import matplotlib.pyplot as plt
import argparse
sys.path.append('../../../../../MCMC/mcmc_mass_age/samples/')
from PercentileClass import PercentileAge


parser = argparse.ArgumentParser()

parser.add_argument('-s',
                    action='store_const',
                    dest='save',
                    const='save',
                    help='save plot'
                    )

parser.add_argument('-p',
                    action='store_const',
                    dest='plot',
                    const='plot',
                    help='show plot'
                    )

parser.add_argument('-l',
                    action='store',
                    dest='system',
                    help='select a system'
                    )



args=parser.parse_args()



system = args.system
Age=PercentileAge(system)


percentiles=['1','2','3','4','5','10','20','30','40','50']
for per in percentiles:

    age=round(Age(int(per)),3)
    if age<1:
        age=age*1e3
        approx_age=str(age)+'Myr'
    else:approx_age=str(age)+'Gyr'
    spinvlogQfilename='System_'+system+'/SpinLogQ_'+per+'.txt'
    p=[]
    q=[]
    with open(spinvlogQfilename,'r') as f:
        next(f)
        for lines in f:
            x=lines.split()
            if numpy.logical_and(abs(float(x[6]))<1e-3,
                                 abs(float(x[7]))<1e-3):
                q=numpy.append(q,float(x[0]))
                p=numpy.append(p,float(x[1]))
    data=zip(q,p)
    data=sorted(data,key=lambda tup: tup[0])
    q_values=[]
    p_values=[]
    for t in data:
        q_values=numpy.append(q_values,t[0])
        p_values=numpy.append(p_values,t[1])
    qmin=min(q_values)
    qmax=max(q_values)
    q_array=numpy.linspace(qmin,qmax,10000)
    p_interpolated=numpy.interp(q_array,q_values,p_values)
    plt.plot(q_array,p_interpolated,label='age='+approx_age)

with open('../../SpinlogQCatalog_el0.4.txt','r') as f:
    for i,lines in enumerate(f):
        if i==0:
            x=lines.split()
            spin_index=x.index('Pspin')
            orbital_index=x.index('Porb')
        else:
            x=lines.split()
            if x[0]==system:
                Pspin=float(x[spin_index])
                Porb=float(x[orbital_index])
                break

plt.hlines(Pspin,5,10,linestyles='dashed',label='Pspin')
plt.hlines(Porb,5,10,linestyles='dotted',label='Porb')
f='xx-small'
plt.legend(loc='upper left',fontsize=f)
plt.title('Spin vs logQ for System '+system)
plt.xlabel('logQ')
plt.ylabel('Spin (days)')
plt.savefig('SpinVsLogQ_'+system+'.pdf')
if args.plot:plt.show()
if args.save:plt.savefig('SpinVsLogQ_'+system+'.pdf')
