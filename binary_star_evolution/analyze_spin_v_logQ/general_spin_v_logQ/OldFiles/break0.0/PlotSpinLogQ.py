import matplotlib.pyplot as plt
import numpy
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('index',help='select system to run')
args = parser.parse_args()

system=args.index


Filename='SpinLogQ_'+system+'_test.txt'

logQ=[]
spin=[]

with open(Filename,'r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        try:
            logQ.append(float(x[0]))
            spin.append(float(x[1]))
        except:
            continue

plt.scatter(logQ,spin)
plt.axhline(y=14.077,color='r',label='Observed Primary Star Spin Period')
plt.axhline(y=15.927,color='k',label='Observed Orbital Period')
plt.xlabel('logQ')
plt.ylabel('Spin Period (days)')
plt.legend()
plt.show()
