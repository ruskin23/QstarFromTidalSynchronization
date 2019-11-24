import matplotlib.pyplot as plt
from numpy import log10
import sys

for k,c in zip([0.5,1.0,1.5,2.0,2.5,3.0],['r','g','b','c','m','y']):
    q=[]
    s=[]
    with open('break'+str(k)+'/SpinLogQ_WithBreaks_'+sys.argv[1]+'.txt','r') as f:
        next(f)
        for lines in f:

            x=lines.split()
            if float(x[7])<1e-5:
                q.append(float(x[0]))
                s.append(log10(float(x[1])))
    plt.plot(q,s,color=c)

plt.axhline(y=log10(12.473),label='PSpinCurrent',color='k')
plt.axhline(y=log10(10.266),label='PorbCurrent',color='darkblue')
plt.legend()
plt.show()


