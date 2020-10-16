import numpy
import matplotlib.pyplot as plt


filename='../Adaptive/gtt/AcceptedParameters.txt'

e=[]
i=[]
with open(filename,'r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        i=numpy.append(i,int(x[0]))
        e=numpy.append(e,float(x[2]))


plt.scatter(i,e)
plt.show()
