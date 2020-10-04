import matplotlib.pyplot as plt
import numpy


p=[]
e=[]
spin=[]
r=[]
with open('SpinlogQCatalog_el0.4.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        if float(x[6])/float(x[12])<10:
            p.append(float(x[6]))
            e.append(float(x[8]))
            spin.append(float(x[12]))
            r.append(float(x[6])/float(x[12]))


plt.axhline(y=1, color='r', linestyle='-')
plt.scatter(p,r)
plt.xscale('log')
plt.xlabel('Orbital Period (days)')
#plt.ylabel('eccentricity')
plt.ylabel('(Orbital Period)/(Spin Period)')
#plt.show()
plt.savefig('PvRCatalog.pdf')


with open('SpinlogQCatalog_el0.4.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        system_main=x[0]
        system_sol=[]
        with open('SolutionFileBreaks0.0.txt','r') as fy:
            next(fy)
            for ylines in fy:
                y=ylines.split()
                system_sol.append(y[0])
        if system_main in system_sol:print(system_main)

