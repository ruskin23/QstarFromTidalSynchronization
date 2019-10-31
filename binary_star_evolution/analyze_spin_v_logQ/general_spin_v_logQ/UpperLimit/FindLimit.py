import os

directory=os.getcwd()
print(directory)
systems=[]
for files in os.listdir(directory):
    x=files.split('_')
    if x[0]=='SpinLogQ':
        systems.append(x[1])
print(systems)

with open('spin_vs_logQ_systems_0.2.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        for s in systems:
            if s==x[0]:
                spinfile='SpinLogQ_'+s+'_test.txt'
                with open(spinfile,'r') as f1:
                    next(f1)
                    for lines1 in f1:
                        y=lines1.split()
                    if float(y[7])<1e-3:print(s)
