import numpy
s=[]
m=[]
with open('SpinlogQCatalog_el0.4.txt','r') as f:
    next(f)
    for lines in f :
        x=lines.split()
        s.append(int(x[0]))

with open('SolutionFileBreaks0.0.txt','r') as fm:
    next(fm)
    for lines in fm:
        y=lines.split()
        m.append(int(y[0]))

set1=set(s)
set2=set(m)
missing=list(sorted(set1-set2))
kmissing=[str(i) for i in missing]
print('Missing Systems:')
print('number of systems= ',len(missing))
print(kmissing)
print('system numbers= ',' '.join(kmissing))


m.sort()
print('System with solutions:')
print('number of systems= ',len(m))
km=[str(i) for i in m]
print('system numbers= ',','.join(km))


s.sort()
print('Total systems:')
print('number of systems= ',len(s))
k=[str(i) for i in s]
print('system numbers= ',' '.join(k))
