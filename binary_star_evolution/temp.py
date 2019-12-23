import numpy
import matplotlib.pyplot as plt


omegaRef=6.28318531e-46
values=[]
with open('temp.txt','r') as f:
    for lines in f:
        x=lines.split(',')

for v in x:
    if v!=',':
        values.append(float(v))

values=numpy.array(values)
#ratio=numpy.where(values<=omegaRef,1,values)

age=numpy.linspace(0,2.4,len(values))

plt.semilogx(age,values)
plt.show()
