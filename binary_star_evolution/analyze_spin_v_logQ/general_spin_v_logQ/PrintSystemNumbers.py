s=[]

with open('SpinlogQCatalog_el0.4.txt','r') as f:
    next(f)
    for lines in f :
        x=lines.split()
        s.append(x[0])

print(' '.join(s))
