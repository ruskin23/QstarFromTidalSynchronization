import os

#directory=os.getcwd()
directory='/home/ruskin/projects/QstarFromTidalSynchronization/binary_star_evolution/analyze_spin_v_logQ/general_spin_v_logQ/break2.0/'
print(directory)
systems=[]
for files in os.listdir(directory):
    x=files.split('_')
    if x[0]=='SpinLogQ':
        y=x[2].split('.')
        systems.append(y[0])
print(systems)

with open('spin_vs_logQ_systems_0.2.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        for s in systems:
            if s==x[0]:
                #spinfile='SpinLogQ_'+s+'_test.txt'
                spinfile=directory+'SpinLogQ_WithBreaks_'+s+'.txt'
                with open(spinfile,'r') as f1:
                    next(f1)
                    for lines1 in f1:
                        y=lines1.split()
                        break
                    if float(y[5])<1e-3:print(s)
