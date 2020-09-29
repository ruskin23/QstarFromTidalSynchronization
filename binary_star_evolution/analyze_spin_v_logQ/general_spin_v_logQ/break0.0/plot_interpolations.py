import sys
import numpy
import matplotlib.pyplot as plt
import argparse


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
spinvlogQfilename='SpinLogQ_'+system+'.txt'
p=[]
q=[]
with open(spinvlogQfilename,'r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        if x[0]!='logQ':
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
    plt.plot(q_array,p_interpolated)

with open('../SpinlogQCatalog_el0.4.txt','r') as f:
    for i,lines in enumerate(f):
        if i==0:
            x=lines.split()
            spin_index=x.index('Pspin')
            orbital_index=x.index('Porb')
        else:
            x=lines.split()
            if x[0]==system:
                kic=x[1]
                Pspin=float(x[spin_index])
                Porb=float(x[orbital_index])
                break


plt.hlines(Pspin,5,12,linestyles='dashed',label='Pspin')
plt.hlines(Porb,5,12,linestyles='dotted',label='Porb')
plt.legend(loc='upper left')
plt.title('Spin vs logQ for KIC'+kic)
plt.xlabel('logQ')
plt.ylabel('Spin (days)')
if args.plot:plt.show()
if args.save:plt.savefig('SpinVsLogQ_'+system+'.eps')
