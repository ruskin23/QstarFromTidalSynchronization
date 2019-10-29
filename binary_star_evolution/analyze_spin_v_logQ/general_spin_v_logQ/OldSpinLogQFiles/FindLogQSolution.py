import numpy

with open('CheckLogQSolution.txt','r') as f:
    for line in f:
        system_array=line.split()
print(system_array)

with open('FindLogQSolution.txt','w') as f:
    for system in system_array:
        print('System = ',system)
        with open('spin_vs_logQ_systems_0.2.txt','r') as f2:
            next(f2)
            for lines in f2:
                x=lines.split()
                if x[0]==system:
                    Porb=float(x[6])
                    Pspin=float(x[12])
                    break
        print('Pspin = ', Pspin)
        print('Porb = ', Porb)

        SpinFilename='SpinLogQ_'+system+'.txt'
        print('SpinFilename = ', SpinFilename)
        with open(SpinFilename,'r') as f1:
            next(f1)
            q=[]
            p=[]
            for lines in f1:
                x=lines.split()
                q.append(float(x[0]))
                p.append(float(x[1]))

            q=numpy.array(q)
            p=numpy.array(p)

            print('LogQArray = ',q)
            print('PspinArray = ',p)


            LogQArray=numpy.linspace(q[0],q[-1],1000)
            PspinInterpolated=numpy.interp(LogQArray,q,p)
            PspinDifference=PspinInterpolated-Pspin
            print(LogQArray)
            print(PspinInterpolated)
            print(PspinDifference)
            zero_crossing = numpy.where(numpy.diff(numpy.sign(PspinDifference)))[0]

            print('Zero_corssing = ',zero_crossing)

            PspinSol=PspinInterpolated[zero_crossing][0]

            print('Solution = ',PspinSol)
            spin_frequency=2*numpy.pi/PspinSol
            orbital_frequency=2*numpy.pi/Porb
            tidal_frequency=2*(spin_frequency-orbital_frequency)
            tidal_period=2*numpy.pi/(tidal_frequency)

            f.write(system+'\t'+
                    repr(LogQArray[zero_crossing][0])+'\t'+
                    repr(PspinSol)+'\t'+
                    repr(Porb)+'\t'+
                    repr(tidal_frequency)+'\t'+
                    repr(tidal_period)+'\n')







