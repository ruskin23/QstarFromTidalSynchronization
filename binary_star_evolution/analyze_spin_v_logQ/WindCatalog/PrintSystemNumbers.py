import numpy
s=[]
with open('NewCatalog.txt','r') as f:
    next(f)
    for lines in f :
        x=lines.split()
        s.append(int(x[0]))

s.sort()
print('Total systems:')
print('number of systems= ',len(s))
k=[str(i) for i in s]
print('system numbers= ',' '.join(k))
