import numpy

logQ =[]
Pspin= []
with open('spin_vs_logQ_KIC8543278.txt','r') as f:
    next(f)
    next(f)
    next(f)
    for line in f:
        data = line.split()
        logQ.append(float(data[0]))
        Pspin.append(float(data[1]))

#logQ=numpy.asarray(logQ)
#Pspin=numpy.asarray(Pspin)


print(type(logQ))

logQ_array = numpy.linspace(7.0,8.0,50)
logQ_interp = numpy.interp(logQ_array,logQ,Pspin)


for i in range(len(logQ_array)-1):
    print(logQ_array[i])
    print(logQ_interp[i])
