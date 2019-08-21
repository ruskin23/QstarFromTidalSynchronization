import numpy
import os

directory='/home/ruskin/projects/QstarFromTidalSynchronization/binary_star_evolution/analyze_spin_v_logQ/general_spin_v_logQ/first_mass_sol/e_range_0.2_0.4/test_format'

count = 0
total = 0

fname = []

for filename in os.listdir(directory):
    if filename.endswith('.txt'):
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
        spin_diff = spin-observed_spin
        if spin_diff[0]*spin_diff[5]<0:
            fname.append(filename)
            count=count+1




for files in fname:
    logQ=[]
    Pspin=[]
    with open(files,'r') as f:
        for i,line in enumerate(f):
            if i==2:
                x=line.split()
                observed_spin=float(x[6])
            if i>2:
                data = line.split()
                logQ.append(float(data[0]))
                Pspin.append(float(data[1]))


    logQ_array = numpy.linspace(5.0,10.0,10000)
    Pspin_interp = numpy.interp(logQ_array,logQ,Pspin)
    Pspin_diff = Pspin_interp-observed_spin

    zero_crossing = numpy.where(numpy.diff(numpy.sign(Pspin_diff)))[0]
    for i in zero_crossing:
        print(logQ_array[i])
        print(Pspin_interp[i])


