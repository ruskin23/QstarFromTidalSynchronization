from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt


pi=[]
ei=[]
pf=[]
ef=[]

with open('temp.txt','r') as f:
    for lines in f:
        x=lines.split()
        if x[0]=='Trying':
            pi.append(float(x[3][0:-2]))
            ei.append(float(x[5][1:]))
        if x[0]=='Final':
            if x[1]=='Eccentricity':
                ef.append(float(x[3]))
            if x[1]=='Porb':
                pf.append(float(x[3]))

print(len(ei))
print(len(pi))
print(len(ef))
print(len(pf))

with open('temp2.txt','w') as f:
    f.write('pi'+'\t'+
            'ei'+'\t'+
            'pf'+'\t'+
            'ef'+'\n')
    for i in range(len(pi)):
        f.write(repr(pi[i])+'\t'+
                repr(ei[i])+'\t'+
                repr(pf[i])+'\t'+
                repr(ef[i])+'\n')


print(len(ei))
print(len(pi))
print(len(ef))


fig=plt.figure()
ax=plt.axes(projection='3d')
ax.scatter3D(pi,ei,ef,c=ef)
plt.show()
