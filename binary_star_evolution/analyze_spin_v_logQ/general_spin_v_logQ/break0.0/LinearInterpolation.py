import numpy
import matplotlib.pyplot as plt
import sys

system=sys.argv[1]
spinfile='PspinInterpolated/SpinlogQ_'+system+'.txt'

with open(spinfile,'r') as f:
        next(f)
        p=[]
        q=[]
        for lines in f:
            x=lines.split()
            q.append(float(x[0]))
            p.append(float(x[1]))

        q=numpy.array(q)
        p=numpy.array(p)

        LogQArray=numpy.linspace(q[0],q[-1],1000)
        PspinInterpolated=numpy.interp(LogQArray,q,p)

plt.scatter(q,p)
plt.plot(LogQArray,PspinInterpolated)
plt.show()
