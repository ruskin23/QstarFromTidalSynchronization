import matplotlib.pyplot as plt
import numpy

p=[]
e=[]
d=[]
pmass=[]
smass=[]
age=[]
feh=[]
system=[]

with open('windemuth_orbital_raw.txt','r') as f:
    next(f)
    next(f)
    for i,lines in enumerate(f):
        x=lines.split()
        if float(x[1])<30:
            system=numpy.append(system,i)
            with open('windemuth_stellar_raw.txt','r') as f2:
                next(f2)
                next(f2)
                for lines2 in f2:
                    y=lines2.split()
                    if y[0]==x[0]:
                        feh=numpy.append(feh,float(y[1]))
                        pmass=numpy.append(pmass,float(y[7]))
                        smass=numpy.append(smass,float(y[10]))
                        age=numpy.append(age,float(y[4]))
                        break
            p=numpy.append(p,float(x[1]))
            e=numpy.append(e,numpy.sqrt(float(x[7])**2 + float(x[10])**2))
            d=numpy.append(d,numpy.log10(float(x[16])))


plt.scatter(p,e)
plt.xscale('log')
plt.yscale('log')
plt.show()


