import numpy
import matplotlib.pyplot as plt

pi=[]
pf=[]

with open('temp.txt','r') as f:
    for lines in f:
        x=lines.split()
        print(x)
        pi.append(numpy.log10(float(x[0])))
        pf.append(numpy.log10(float(x[1])))
        #if x[0]=='Trying':
        #    pi.append(float(x[3][0:-2]))
        #if x[0]=='Final':
        #    pf.append(float(x[3]))


plt.scatter(pi,pf)
plt.show()
