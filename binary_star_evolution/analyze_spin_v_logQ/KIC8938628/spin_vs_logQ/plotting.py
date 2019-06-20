import matplotlib.pyplot as plt
import numpy

logq=[]
spin=[]

with open('check_sol_1.0.txt','r') as f:
    for line in f:
        x=line.split()
        logq.append(float(x[0]))
        spin.append(float(x[1]))


logq=numpy.array(logq)
spin=numpy.array(spin)

spin_max = 19.304 + 2.689
spin_min = 19.304 - 2.689

plt.scatter(logq,spin)
plt.axhline(y=spin_max, color='g', linestyle='--',label='+1$\sigma$')
plt.axhline(y=spin_min, color='g', linestyle='--',label='-1$\sigma$')
plt.xlabel('logQ', fontsize=18)
plt.legend()
plt.ylabel('spin period (days)', fontsize=16)
plt.show()
