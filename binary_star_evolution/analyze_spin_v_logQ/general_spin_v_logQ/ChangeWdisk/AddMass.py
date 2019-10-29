import numpy


system_array=[]
with open('FindLogQSolution.txt','r') as f:
    for lines in f:
        x=lines.split()
        system_array.append(x[0])


mass_array=[]
with open('spin_vs_logQ_systems_0.2.txt','r') as f:
    for system in system_array:
        for lines in f:
            x=lines.split()
            at_sys=x[0]
            if at_sys==system:
                mass_array.append(x[15])
                break



with open('FindLogQSolution_withMass.txt','w') as f:
    with open('FindLogQSolution.txt','r') as f1:
        for i,lines in enumerate(f1):
            x=lines.split()
            x.append(mass_array[i])
            y='\t'.join(x)
            f.write(y+'\n')



