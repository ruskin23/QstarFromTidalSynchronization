import matplotlib.pyplot as plt
import numpy

q=[]
s=[]
with open('SpinLogQ_1_test.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        q.append(float(x[0]))
        s.append(float(x[1]))

Pspin=10.778
Porb=8.46

plt.axhline(y=Pspin,label='PSpinCurrent',color='r',linestyle='--')
plt.axhline(y=Porb,label='PorbCurrent',color='k',linestyle='--')
plt.plot(q,s)
plt.legend()
plt.xlabel('logQ')
plt.ylabel('Ppsin (days)')
plt.show()
