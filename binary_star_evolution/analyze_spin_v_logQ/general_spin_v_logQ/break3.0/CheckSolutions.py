import os
import numpy

with open('FindingSolution.txt','w') as f:
    f.write('system'+'\t'+
            'logq'+'\t'+
            'pspin'+'\t'+
            'tidal_frequency'+'\t'+
            'tidal_period'+'\n')

def FindSolution(system,
                 spinfile,
                 Pspin,
                 Porb
                 ):

    with open(spinfile,'r') as f:
        next(f)
        p=[]
        q=[]
        for lines in f:
            x=lines.split()
            q.append(float(x[0]))
            p.append(float(x[1]))

        q=numpy.array(q)
        p=numpy.array(p)

        LogQArray=numpy.linspace(q[0],q[-1],1000)
        PspinInterpolated=numpy.interp(LogQArray,q,p)
        PspinDifference=PspinInterpolated-Pspin

        zero_crossing = numpy.where(numpy.diff(numpy.sign(PspinDifference)))[0]

        PspinSol=PspinInterpolated[zero_crossing][0]

        spin_frequency=2*numpy.pi/PspinSol
        orbital_frequency=2*numpy.pi/Porb
        tidal_frequency=2*(spin_frequency-orbital_frequency)
        tidal_period=2*numpy.pi/(tidal_frequency)

        with open('FindingSolution.txt','a') as f1:
            f1.write(system+'\t'+
                    repr(LogQArray[zero_crossing][0])+'\t'+
                    repr(PspinSol)+'\t'+
                    repr(Porb)+'\t'+
                    repr(tidal_frequency)+'\t'+
                    repr(tidal_period)+'\n')


directory=os.getcwd()
print(directory)
systems=[]
for files in os.listdir(directory):
    x=files.split('_')
    if x[0]=='SpinLogQ':
        name=x[2].split('.')
        systems.append(name[0])
print(systems)
with open('spin_vs_logQ_systems_0.2.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        for s in systems:
            if x[0]==s:
                p=[]
                spin=float(x[12])
                porb=float(x[6])
                with open('SpinLogQ_WithBreaks_'+s+'.txt','r') as f1:
                    next(f1)
                    for lines1 in f1:
                        y=lines1.split()
                        p.append(float(y[1])-spin)
                    print('For system: ',s)
                    print(p)
                    p=numpy.array(p)
                    zero_crossing=numpy.where(numpy.diff(numpy.sign(p)))[0]
                    if zero_crossing.size==1:
                        FindSolution(s,
                                     'SpinLogQ_WithBreaks_'+s+'.txt',
                                     spin,
                                     porb
                                     )
                        print(s)

                    break
