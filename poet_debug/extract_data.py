#Exctract data from stopping condition lines.
import matplotlib.pyplot as plt
import os
import numpy
import argparse

parser=argparse.ArgumentParser()
parser.add_argument('condition_no')
args = parser.parse_args()



fname_in='stopping_condition_output_3.txt'
fpath='/mnt/md0/ruskin/QstarFromTidalSynchronization/poet_debug'


in_file=os.path.join(fpath,fname_in)



terms_age = ['Age:','Age:|']
age=[]

terms_condition =['Condition['+args.condition_no+']:','Condition['+args.condition_no+']:|']
condition=[]

with open(in_file,'r') as f:
    for line in f:
        d=line.split()
        if d[0] in terms_age:
            print(line)
            del d[0]
            for item in d:
                if '|' in item:item=item.replace('|','')
                age.append(float(item))
        if len(d)>1:
            a=d[0]+d[1]
            if a in terms_condition:
                for i in range(4):
                    del d[0]
                for item in d:
                    if '|' in item:item=item.replace('|','')
                    condition.append(float(item))

print(len(age))
print(len(condition))

fname_out='tabulated_stopping_condition_'+args.condition_no+'.txt'
out_file=os.path.join(fpath,fname_out)

with open(out_file,'w') as f:
    for a,b in zip(age,condition):
        f.write(repr(a)+'\t'+repr(b)+'\n')

plt.plot(age,condition)
plt.show()
