import sys


system=sys.argv[1]

with open('../SpinlogQCatalog_el0.4.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        if x[0]==system:
            Pspin=float(x[12])
            PspinE=float(x[13])
            break

print('{} {}'.format(Pspin,PspinE))
