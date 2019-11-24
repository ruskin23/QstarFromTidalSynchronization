import os
import numpy

with open('SolutionFile.txt','w') as f:
    f.write('system'+'\t'+
            'logq'+'\t'+
            'pspin'+'\t'+
            'porb'+'\t'+
            'tidal_frequency'+'\t'+
            'tidal_period'+'\t'+
            'logtidal_frequency'+'\t'+
            'logtidal_period'+'\t'+
            'primary_mass'+'\t'+
            'synchronized'+'\t'+
            'circularised'+'\n')

def FindUpperLimit(system,
                   spinfile,
                   spin,
                   porb,
                   primary_mass
                   ):

    circularized=None
    synchronized=None
    directoy='/home/ruskin/projects/QstarFromTidalSynchronization/binary_star_evolution/analyze_spin_v_logQ/general_spin_v_logQ/UpperLimit/NoBreaks/'
    filename=directoy+spinfile

    if os.path.exists(filename):
        with open(directoy+spinfile,'r') as f:
            next(f)
            for lines in f:
                x=lines.split()
            if float(x[7])<1e-3:
                circularized=True
                logq=float(x[0])
                spin_frequency=2*numpy.pi/spin
                orbital_frequency=2*numpy.pi/porb
                tidal_frequency=2*(spin_frequency-orbital_frequency)
                tidal_period=2*numpy.pi/(tidal_frequency)
                with open('SolutionFile.txt','a') as f1:
                    f1.write(system+'\t'+
                             repr(logq)+'\t'+
                             repr(spin)+'\t'+
                             repr(porb)+'\t'+
                             repr(tidal_frequency)+'\t'+
                             repr(tidal_period)+'\t'+
                             repr(numpy.log10(tidal_frequency))+'\t'+
                             repr(numpy.log10(tidal_period))+'\t'+
                             repr(primary_mass)+'\t'+
                             str(synchronized)+'\t'+
                             str(circularized)+'\n')


def FindSolution(system,
                 spinfile,
                 Pspin,
                 Porb,
                 Pspin_error,
                 primary_mass
                 ):

    circularized=None
    synchronized=None

    with open(spinfile,'r') as f:
        next(f)
        p=[]
        q=[]
        for lines in f:
            x=lines.split()
            try:
                q.append(float(x[0]))
                p.append(float(x[1]))
            except:continue

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

        PspinSigmaMin=Pspin-Pspin_error
        PspinSigmaMax=Pspin+Pspin_error

        if Porb>PspinSigmaMin and Porb<PspinSigmaMax:
            synchronized=True
        else: synchronized=None

        with open('SolutionFile.txt','a') as f1:
            f1.write(system+'\t'+
                    repr(LogQArray[zero_crossing][0])+'\t'+
                    repr(PspinSol)+'\t'+
                    repr(Porb)+'\t'+
                    repr(tidal_frequency)+'\t'+
                    repr(tidal_period)+'\t'+
                    repr(numpy.log10(abs(tidal_frequency)))+'\t'+
                    repr(numpy.log10(abs(tidal_period)))+'\t'+
                    repr(primary_mass)+'\t'+
                    str(synchronized)+'\t'+
                    str(circularized)+'\n')


directory=os.getcwd()
print(directory)
systems=[]
for files in os.listdir(directory):
    x=files.split('_')
    if x[0]=='SpinLogQ':
        y=x[2].split('.')
        systems.append(y[0])
print(systems)
with open('spin_vs_logQ_systems.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        for s in systems:
            if x[0]==s:
                p=[]
                spin=float(x[12])
                porb=float(x[6])
                spin_error=float(x[13])
                primary_mass=float(x[14])
                with open('SpinLogQ_WithBreaks_'+s+'.txt','r') as f1:
                    next(f1)
                    lowspin_check=0
                    incomplete_check=0
                    for lines1 in f1:
                        y=lines1.split()
                        try:
                            if y[1]!='spin':
                                p.append(float(y[1])-spin)
                        except:continue
                    if len(p)<6:
                        incomplete_check=1
                        print('Incomplete for system: ',s)
                    p=numpy.array(p)
                    zero_crossing=numpy.where(numpy.diff(numpy.sign(p)))[0]

                    circularized=None
                    synchronized=None

                    if zero_crossing.size==1:
                        FindSolution(s,
                                     'SpinLogQ_WithBreaks_'+s+'.txt',
                                     spin,
                                     porb,
                                     spin_error,
                                     primary_mass
                                     )
                        print('systems with single solution: ',s)
                    if zero_crossing.size==0 and spin/porb>1 and incomplete_check==0 and y[0]>y[-1]:
                        print('systems cannot produce high pspin values: ',s)
                        lowspin_check=1
                    if zero_crossing.size==0 and spin/porb<1 and lowspin_check==0 and incomplete_check==0:
                        FindUpperLimit(s,
                                       'SpinLogQ_WithBreaks_'+s+'.txt',
                                       spin,
                                       porb,
                                       primary_mass)
                        print('system with no solutions: ',s)
                    if zero_crossing.size>1:
                        print('system with mulitple solution: ',s)
                    break
