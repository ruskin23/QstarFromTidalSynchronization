import matplotlib.pyplot as plt
import numpy

porb=[]
pspin=[]
r=[]
e=[]
with open('catalog_KIC.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        porb.append(float(x[5]))
        pspin.append(float(x[11]))
        if float(x[5])/float(x[11])<6:
            r.append(float(x[5])/float(x[11]))
            e.append(float(x[7]))
        if float(x[5])/float(x[11])>1.2 and float(x[7])>0.1 and float(x[7])<0.3:
            print(lines)


e=numpy.array(e)
porb=numpy.array(porb)
pspin=numpy.array(pspin)
r=numpy.array(r)

plt.scatter(e,r)
plt.axhline(y=1)
plt.show()


