import numpy
s=[]

with open('SpinlogQCatalog_el0.4.txt','r') as f:
    next(f)
    for lines in f :
        x=lines.split()
        s.append(int(x[0]))

s.sort()
print(s)
k=[str(i) for i in s]
print(' '.join(k))
