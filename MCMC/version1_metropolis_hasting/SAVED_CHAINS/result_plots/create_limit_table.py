import sys
import matplotlib.pyplot as plt
import numpy
import itertools

system=sys.argv[1]
qdata=[]
adata=[]
with open('AcceptedParameters.txt','r') as f:
    for lines in f:
        x=lines.split()
        try:
            q=float(x[4])
            a=float(x[6])
            qdata=numpy.append(qdata,q)
            adata=numpy.append(adata,a)
        except:continue


qvalues=[(x, len(list(y))) for x, y in itertools.groupby(qdata)]
qvalues=sorted(qvalues, key=lambda tup: tup[0])

value=[]
multiplicity=[]
m=0
for a in qvalues:
    value=numpy.append(value,a[0])
    m=m+a[1]
    multiplicity=numpy.append(multiplicity,m)

multiplicity=multiplicity/max(multiplicity)


percentiles=[0.1,0.33,0.5,0.66,0.9]
value_array=[]
for i in range(5):
    error=1e-10
    while True:
        for v,m in zip(value,multiplicity):
            if abs(m-percentiles[i])<error:
                value_array=numpy.append(value_array,v)
                break
        if len(value_array)<i+1:
            error=error*10
            continue
        else:break

low_percentile=value_array[1]
mean=value_array[2]
up_percentile=value_array[3]

low_limit=value_array[0]
up_limit=value_array[4]

with open('../SpinlogQCatalog_el0.4.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        if x[0]==system:
            porb=float(x[6])
            pspin=float(x[12])
            pspin_error=float(x[13])
            break

pspin_up=pspin+2*pspin_error
pspin_down=pspin-2*pspin_error

if numpy.logical_and(porb>pspin_down,porb<pspin_up):
    tidal_frequency=0
    tidal_period=porb
else:
    tidal_frequency=2*abs((2*numpy.pi/porb)-(2*numpy.pi/pspin))
    tidal_period=2*numpy.pi/tidal_frequency


with open('limits.txt','a') as f:
    print('\nSysten = ',system)
    f.write(system+'\t'+
            repr(mean)+'\t'+
            repr(low_limit)+'\t'+
            repr(up_limit)+'\t'+
            repr(low_percentile)+'\t'+
            repr(up_percentile)+'\t'+
            repr(tidal_frequency)+'\t'+
            repr(tidal_period)+'\t'+
            repr(porb)+'\t'+
            repr(pspin)+'\n')

