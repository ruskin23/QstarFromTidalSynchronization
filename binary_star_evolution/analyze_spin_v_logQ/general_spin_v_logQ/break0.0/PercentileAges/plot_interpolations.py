import sys
import numpy
import matplotlib.pyplot as plt

system = sys.argv[1]
#extra_age=['3.2','3.4','3.6','3.8']
age=['1','2','3','4','5','10','20','30','40','50']
#age=age+extra_age
for a in age:
    spinvlogQfilename='System_'+system+'/SpinLogQ_'+a+'.txt'
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
        print('For age {} q0={} qn={}'.format(a,qmin,qmax))
        q_array=numpy.linspace(qmin,qmax,10000)
        p_interpolated=numpy.interp(q_array,q_values,p_values)
        plt.plot(q_array,p_interpolated,label='age='+a)

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
plt.legend(loc='upper left')
plt.show()
