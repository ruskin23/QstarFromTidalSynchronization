import numpy

system_array=[]
PspinValue=[]
with open('spin_vs_logQ_systems_0.2.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        system_array.append(x[0])
        PspinValue.append(float(x[12]))

with open('CheckLogQSolution.txt','w') as f:
    for spin,system in zip(PspinValue,system_array):
        SpinFilename='SpinLogQ_'+system+'.txt'
        with open(SpinFilename,'r') as f1:
            next(f1)
            p=[]
            for i,lines in enumerate(f1):
                x=lines.split()
                p.append(float(x[1])-spin)
            p=numpy.array(p)
            zero_crossing=numpy.where(numpy.diff(numpy.sign(p)))[0]
            if zero_crossing.size==1:

                f.write(system+'\t')
            elif zero_crossing.size==0:
                print(system)
