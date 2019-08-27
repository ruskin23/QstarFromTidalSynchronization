import numpy
import os

directory='/home/ruskin/projects/QstarFromTidalSynchronization/binary_star_evolution/analyze_spin_v_logQ/general_spin_v_logQ/first_mass_sol/linear_interp/'

count = 0
total = 0

fname = []
fname_no_sol = []
for filename in os.listdir(directory):
    if filename.endswith('.txt'):
        print(filename)
        total=total+1
        spin=[]
        with open(filename,'r') as f:
            for i,line in enumerate(f):
                if i==2:
                    x=line.split()
                    observed_spin=float(x[6])
                if i>2:
                    x=line.split()
                    spin.append(float(x[1]))

        spin = numpy.array(spin)
        print(spin)
        spin_diff = spin-observed_spin
        print(spin_diff)
        if spin_diff[0]*spin_diff[-1]<0:
            fname.append(filename)
            count=count+1
        else:
            fname_no_sol.append(filename)
fm = open('no_sol_files.csv','w')
fm.close()
for x in fname_no_sol:
    with open('no_sol_files.csv','a') as fn:
        fn.write(x +'\n')



logQ_sol=[]
Pspin_sol=[]
filen=[]
Porb=[]

for files in fname:
    logQ=[]
    Pspin=[]
    with open(files,'r') as f:
        for i,line in enumerate(f):
            if i==2:
                x=line.split()
                observed_spin=float(x[6])
                observed_Porb=float(x[5])
                print(observed_Porb)
            if i>2:
                data = line.split()
                logQ.append(float(data[0]))
                Pspin.append(float(data[1]))



    logQ_array = numpy.linspace(logQ[0],logQ[-1],10000)
    Pspin_interp = numpy.interp(logQ_array,logQ,Pspin)
    Pspin_diff = Pspin_interp-observed_spin

    zero_crossing = numpy.where(numpy.diff(numpy.sign(Pspin_diff)))[0]

    for i in zero_crossing:
        logQ_sol.append(logQ_array[i])
        Pspin_sol.append(Pspin_interp[i])
        filen.append(files)
        Porb.append(float(observed_Porb))

print(Porb)


with open('sol_file.csv','w') as f:
    for i in range(len(logQ_sol)):
        print(i)
        spin_frequency=2*numpy.pi/Pspin_sol[i]
        orbital_frequency=2*numpy.pi/Porb[i]
        tidal_frequency=2*(spin_frequency-orbital_frequency)
        f.write(filen[i] + '\t' + repr(logQ_sol[i]) +'\t'+ repr(Pspin_sol[i])
                +'\t'+ repr(Porb[i]) +'\t'+ repr(2*numpy.pi/tidal_frequency) + '\n')
