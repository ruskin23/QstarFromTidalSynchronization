from matplotlib import pyplot as plt
import numpy
s=[]
p=[]
q=[]
with open('FindingSolution.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        s.append(numpy.log10(abs(float(x[4]))))
        p.append(numpy.log10(abs(float(x[5]))))
        q.append(float(x[1]))

with open('NewUpper.txt','w') as f:
    for i in range(len(s)):
        f.write(repr(q[i])+'\t'+repr(s[i])+'\t'+repr(p[i])+'\n')

