from matplotlib import pyplot as plt
import numpy


synchronised={'logQ': [],
              'logfrequency':[],
              'logperiod':[],
              'mass':[]}

circularised={'logQ': [],
              'logfrequency':[],
              'logperiod':[],
              'mass':[]}
exact={'logQ': [],
       'logfrequency':[],
       'logperiod':[],
       'mass':[]}

with open('SolutionFile.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        sync=x[9]
        circ=x[10]
        if sync=='True' and circ=='None':
            synchronised['logQ'].append(float(x[1]))
            synchronised['logfrequency'].append(float(x[6]))
            synchronised['logperiod'].append(float(x[7]))
            synchronised['mass'].append(float(x[8]))

        if sync=='None' and circ=='True':
            circularised['logQ'].append(float(x[1]))
            circularised['logfrequency'].append(float(x[6]))
            circularised['logperiod'].append(float(x[7]))
            circularised['mass'].append(float(x[8]))

        if sync=='None' and circ=='None':
            exact['logQ'].append(float(x[1]))
            exact['logfrequency'].append(float(x[6]))
            exact['logperiod'].append(float(x[7]))
            exact['mass'].append(float(x[8]))



#plt.scatter(synchronised['logfrequency'],synchronised['logQ'],marker='v',label='synchronised')
#plt.scatter(circularised['logfrequency'],circularised['logQ'],marker='v',label='circularised')
#plt.scatter(exact['logfrequency'],exact['logQ'],marker='o',label='exact')
#plt.legend()
#plt.savefig('TidalFrequencyNoBreaks.eps')


#plt.scatter(synchronised['logperiod'],synchronised['logQ'],marker='v',label='synchronised')
#plt.scatter(circularised['logperiod'],circularised['logQ'],marker='v',label='circularised')
#plt.scatter(exact['logperiod'],exact['logQ'],marker='o',label='exact')
#plt.legend()
#plt.savefig('TidalPeriodNoBreaks.eps')


plt.scatter(synchronised['mass'],synchronised['logQ'],marker='v',label='synchronised')
plt.scatter(circularised['mass'],circularised['logQ'],marker='v',label='circularised')
plt.scatter(exact['mass'],exact['logQ'],marker='o',label='exact')
plt.legend()
plt.savefig('PrimaryMassNoBreaks.eps')
