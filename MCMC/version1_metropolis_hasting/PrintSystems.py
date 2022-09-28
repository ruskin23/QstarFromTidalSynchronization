

s=[]
ns=[]
a=[]
with open('SolutionFileBreaks0.0.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        a.append(x[0])
        if x[10]=='None':
            s.append(x[0])
        if x[10]=='True':
            ns.append(x[0])

print('Non Synchronous Systems: ',len(s))
print(' '.join(s))
print('Synchronous Systems: ',len(ns))
print(' '.join(ns))
print('Total: ',len(ns+s))
print(a)