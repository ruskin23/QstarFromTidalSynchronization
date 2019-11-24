from matplotlib import pyplot as plt
import numpy


synchronised={'logQ': [],
              'logfrequency':[],
              'logperiod':[],
              'mass':[],
              'age':[]}

exact={'logQ': [],
       'logfrequency':[],
       'logperiod':[],
       'mass':[],
       'age':[]}

with open('SolutionFile.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        sync=x[10]
        if sync=='True':
            synchronised['logQ'].append(float(x[1]))
            synchronised['logfrequency'].append(float(x[6]))
            synchronised['logperiod'].append(float(x[7]))
            synchronised['mass'].append(float(x[8]))
            synchronised['age'].append(float(x[9]))

        if sync=='None':
            exact['logQ'].append(float(x[1]))
            exact['logfrequency'].append(float(x[6]))
            exact['logperiod'].append(float(x[7]))
            exact['mass'].append(float(x[8]))
            exact['age'].append(float(x[9]))


plt.scatter(synchronised['logfrequency'],synchronised['logQ'],marker='v',label='synchronised limit')
plt.scatter(exact['logfrequency'],exact['logQ'],marker='o',label='exact')
plt.legend()
plt.xlabel('Tidal Frequency')
plt.ylabel('logQ')
plt.show()

#plt.scatter(synchronised['mass'],synchronised['logQ'],marker='v',label='synchronised')
#plt.scatter(circularised['mass'],circularised['logQ'],marker='v',label='circularised')
#plt.scatter(exact['mass'],exact['logQ'],marker='o',label='exact')
#plt.legend()
#plt.savefig('PrimaryMassNoBreaks.eps')
