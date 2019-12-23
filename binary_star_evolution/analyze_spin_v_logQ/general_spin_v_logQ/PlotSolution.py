from matplotlib import pyplot as plt
import numpy


synchronised={'logQ': [],
              'logfrequency':[],
              'logperiod':[],
              'mass':[],
              'age':[],
              'eccentricity':[]}

exact={'logQ': [],
       'logfrequency':[],
       'logperiod':[],
       'mass':[],
       'age':[],
       'eccentricity':[]}

with open('SpinVslogQ_WithErrorBars.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        sync=x[14]
        if sync=='True':
            synchronised['logQ'].append(float(x[1]))
            synchronised['logfrequency'].append(float(x[11]))
            synchronised['logperiod'].append(abs(float(x[10])))
            synchronised['mass'].append(float(x[12]))
            synchronised['age'].append(float(x[13]))
            synchronised['eccentricity'].append(float(x[6]))
        if sync=='None':
            exact['logQ'].append(float(x[1]))
            exact['logfrequency'].append(float(x[11]))
            exact['logperiod'].append(abs(float(x[10])))
            exact['mass'].append(float(x[12]))
            exact['age'].append(float(x[13]))
            exact['eccentricity'].append(float(x[6]))

plt.semilogx(synchronised['logperiod'],synchronised['logQ'],'v',label='synchronised limit')
plt.semilogx(exact['logperiod'],exact['logQ'],'o',label='exact')
plt.legend()
plt.xlabel('Tidal Frequency')
plt.ylabel('logQ')
plt.show()

#plt.scatter(synchronised['mass'],synchronised['logQ'],marker='v',label='synchronised')
#plt.scatter(circularised['mass'],circularised['logQ'],marker='v',label='circularised')
#plt.scatter(exact['mass'],exact['logQ'],marker='o',label='exact')
#plt.legend()
#plt.savefig('PrimaryMassNoBreaks.eps')
