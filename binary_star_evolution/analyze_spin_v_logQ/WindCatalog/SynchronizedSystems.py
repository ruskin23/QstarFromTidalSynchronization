import numpy

s=[]
ts=[]
with open('NewCatalog.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        ts.append(x[0])
        Porb=float(x[2])
        PspinValue=float(x[14])
        PspinError=float(x[15])
        PspinMax=PspinValue+PspinError
        PspinMin=PspinValue-PspinError
        if numpy.logical_and(Porb>PspinMin,
                             Porb<PspinMax):
            s.append(x[0])

print('Total Systems = ',len(ts))
print(' '.join(ts))
print('Synchronized Systems = ',len(s))
print(' '.join(s))
