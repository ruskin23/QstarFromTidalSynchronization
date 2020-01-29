import os
import numpy
import sys

breakPower=sys.argv[1]

with open('SolutionFileBreaks'+breakPower+'.txt','w') as f:
    f.write('system'+'\t'+
            'logq'+'\t'+
            'pspin'+'\t'+
            'porb'+'\t'+
            'tidal_frequency'+'\t'+
            'tidal_period'+'\t'+
            'logtidal_frequency'+'\t'+
            'logtidal_period'+'\t'+
            'primary_mass'+'\t'+
            'age'+'\t'+
            'synchronized'+'\n')

def FindSynchronizationLimit(key,
                             logQValues,
                             PspinValues,
                             PspinObserved,
                             PspinError,
                             Porb,
                             PrimaryMass,
                             Age):

    synchronized=True

    LogQArray=numpy.linspace(logQValues[0],logQValues[-1],1000)
    PspinInterpolated=numpy.interp(LogQArray,logQValues,PspinValues)

    with open('break'+breakPower+'/PspinInterpolated/SpinlogQ_'+key+'.txt','w') as f:
        for q,s in zip(LogQArray,PspinInterpolated):
            f.write(repr(q)+'\t'+repr(s)+'\t'+repr(s-Porb)+'\n')

    PspinDiff=PspinInterpolated-Porb

    error=1e-10
    i=0

    while error<1e-3:
        error=abs(PspinDiff[i])
        i=i+1

    LimitIndex=i

    LogQLimit=LogQArray[LimitIndex]


    if PspinObserved==Porb:PspinSol=PspinObserved+PspinError
    else:PspinSol=PspinObserved

    spin_frequency=2*numpy.pi/PspinSol
    orbital_frequency=2*numpy.pi/Porb
    tidal_frequency=2*(spin_frequency-orbital_frequency)
    tidal_period=2*numpy.pi/(tidal_frequency)

    with open('SolutionFileBreaks'+breakPower+'.txt','a') as f1:
        f1.write(key+'\t'+
                repr(LogQLimit)+'\t'+
                repr(PspinSol)+'\t'+
                repr(Porb)+'\t'+
                repr(tidal_frequency)+'\t'+
                repr(tidal_period)+'\t'+
                repr(numpy.log10(abs(tidal_frequency)))+'\t'+
                repr(numpy.log10(abs(tidal_period)))+'\t'+
                repr(PrimaryMass)+'\t'+
                repr(Age)+'\t'+
                str(synchronized)+'\n')





def FindSolution(key,
                 logQValues,
                 PspinValues,
                 PspinObserved,
                 PspinError,
                 Porb,
                 PrimaryMass,
                 Age):

    synchronized=None

    LogQArray=numpy.linspace(logQValues[0],logQValues[-1],10000)
    PspinInterpolated=numpy.interp(LogQArray,logQValues,PspinValues)

    with open('break'+breakPower+'/PspinInterpolated/SpinlogQ_'+key+'.txt','w') as f:
        for q,s in zip(LogQArray,PspinInterpolated):
            f.write(repr(q)+'\t'+repr(s)+'\n')

    PspinDiff=PspinInterpolated-PspinObserved

    zero_crossing=numpy.where(numpy.diff(numpy.sign(PspinDiff)))[0]
    if len(zero_crossing)>1:
        print('Multliple Solution for :',key)

    PspinSol=PspinInterpolated[zero_crossing][0]

    PspinMin=PspinObserved-PspinError
    PspinMax=PspinObserved-PspinError

    if numpy.logical_and(Porb>PspinMin,
                         Porb<PspinMax):
        synchronized=True

    spin_frequency=2*numpy.pi/PspinSol
    orbital_frequency=2*numpy.pi/Porb
    tidal_frequency=2*(spin_frequency-orbital_frequency)
    tidal_period=2*numpy.pi/(tidal_frequency)

    with open('SolutionFileBreaks'+breakPower+'.txt','a') as f1:
        f1.write(key+'\t'+
                repr(LogQArray[zero_crossing][0])+'\t'+
                repr(PspinSol)+'\t'+
                repr(Porb)+'\t'+
                repr(tidal_frequency)+'\t'+
                repr(tidal_period)+'\t'+
                repr(numpy.log10(abs(tidal_frequency)))+'\t'+
                repr(numpy.log10(abs(tidal_period)))+'\t'+
                repr(PrimaryMass)+'\t'+
                repr(Age)+'\t'+
                str(synchronized)+'\n')



def SearchSolution(SystemDict):

    for key in SystemDict:

        with open('SpinlogQCatalog_el0.4.txt','r') as f:
            next(f)
            for lines in f:
                x=lines.split()
                at_system=x[0]
                if at_system==key:
                    print('Searchin for system = ',x[0])
                    PspinObserved=float(x[12])
                    Porb=float(x[6])
                    PspinError=float(x[13])
                    PrimaryMass=float(x[15])
                    Age=float(x[16])
                    PspinValues=numpy.array(SystemDict[key]['spin'])
                    logQValues=numpy.array(SystemDict[key]['logQ'])
                    PspinDiff=PspinValues-PspinObserved
                    if PspinDiff[0]*PspinDiff[-1]<0:
                        print('Solution Found for = ', key)
                        FindSolution(key,
                                     logQValues,
                                     PspinValues,
                                     PspinObserved,
                                     PspinError,
                                     Porb,
                                     PrimaryMass,
                                     Age)
                    elif abs(PspinValues[0]-Porb)<0.01:
                        if numpy.logical_and(Porb>PspinObserved-PspinError,
                                             Porb<PspinObserved+PspinError):
                            print('Limit for = ', key)
                            FindSynchronizationLimit(key,
                                                     logQValues,
                                                     PspinValues,
                                                     PspinObserved,
                                                     PspinError,
                                                     Porb,
                                                     PrimaryMass,
                                                     Age)
                        else:
                            print('System Number: ', key)
                    else:
                        print('No Solution for : ', key)
                    break

directory='break'+breakPower+'/'
systems=[]
for files in os.listdir(directory):
    x=files.split('_')
    if x[0]=='SpinLogQ':
        y=x[1].split('.')
        systems.append(y[0])

print(systems)


SystemDict=dict()
for s in systems:
    number=str(s)
    lgq=[]
    spin=[]
    if s=='35':continue
    with open('break'+breakPower+'/SpinLogQ_'+s+'.txt','r') as f:
        if s=='59':print('Found system 59')
        for i,lines in enumerate(f):
            x=lines.split()
            if x[0]=='logQ':continue
            if abs(float(x[7]))<1e-5:
                lgq.append(float(x[0]))
                spin.append(float(x[1]))

    if i<4:
        print('Incomplete for system = ',s)
    else:
        SystemDict[number]=dict()
        SystemDict[number]['logQ']=lgq
        SystemDict[number]['spin']=spin

print(SystemDict['59'])
SearchSolution(SystemDict)
