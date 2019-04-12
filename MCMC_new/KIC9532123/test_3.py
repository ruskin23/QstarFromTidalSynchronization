import matplotlib.pyplot as plt
import numpy


values = ['4.6','4.7','4.8','4.9','5.0','5.1','5.2','5.3','5.4','5.5']

for v in values:
    fname = 'spin_'+v+'.txt'
    e_i=[]
    e_f=[]
    with open(fname,'r') as f:
        next(f)
        for lines in f:
            x=lines.split()
            e_i.append(float(x[0]))
            e_f.append(float(x[1]))
    plt.scatter(e_i,e_f,label=fname,marker='x')
plt.legend(loc='upper left')
plt.show()
