import matplotlib.pyplot as plt
import argparse
import numpy

parser = argparse.ArgumentParser()
parser.add_argument('-c',dest='condition',help='select system to run')
args = parser.parse_args()

cond=2*int(args.condition)

with open('checking_output_134.txt','r') as f:
    for i,lines in enumerate(f):
        if i==0:
            x=lines.split()
            age=[x[k] for k in range(1,len(x))]
            print(len(age))
        elif i==cond:
            x=lines.split()
            condition=[x[k] for k in range(4,len(x))]
            print(len(condition))


plt.plot(age,condition)
plt.show()
