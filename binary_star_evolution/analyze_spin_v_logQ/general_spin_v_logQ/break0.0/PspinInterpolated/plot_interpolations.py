import os
import sys
import re
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

#directory='System_'+system
#regex=re.compile(r'\d+')
#percentiles=[]

#for filename in os.listdir(directory):
#    if filename.endswith(".txt"):
#        p=regex.findall(filename)
#        if len(p)==1:
#            percentiles.append(p[0])
#        elif len(p)==2:
#            p1=p[0]
#            p2=p[1]
#            a=(p1+'.'+p2)
#            print(a)
#            percentiles.append(a)
#        else:continue
#percentiles.sort(key=float)

#percentiles=['1','2','3','4','5','10','20','30','40','50']
#print(percentiles)
spinvlogQfilename='SpinlogQ_'+system+'.txt'
p=[]
q=[]
with open(spinvlogQfilename,'r') as f:
    for lines in f:
        x=lines.split()
        if x[0]=='logQ':continue
        q=numpy.append(q,float(x[0]))
        p=numpy.append(p,float(x[1]))
plt.plot(q,p)

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
if args.plot:plt.show()
if args.save:plt.savefig('SpinVsLogQ_'+system+'.pdf')
