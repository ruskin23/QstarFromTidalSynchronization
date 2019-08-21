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

fname ='last_stopping_condition.txt'

with open(fname,'r') as f:
    for i,line in enumerate(f):
        value = line.split()
        if i==0:age=value
        if i==1:condition0=value
        if i==2:condition1=value
        if i==3:condition2=value
        if i==4:condition3=value
        if i==5:condition4=value
        if i==6:condition5=value
        if i==7:condition6=value

remove_con0=['Condition[', '0]:', '\\', 'ze>']
remove_con1=['Condition[', '1]:', '\\', 'ze>']
remove_con2=['Condition[', '2]:', '\\', 'ze>']
remove_con3=['Condition[', '3]:', '/', 'ze>']
remove_con4=['Condition[', '4]:', '/', 'ze>']
remove_con5=['Condition[', '5]:', '\\', 'ze>']
remove_con6=['Condition[', '6]:', '/', 'ze>']
remove_age=['Age:']

for items in remove_age:
    age.remove(items)
for items in remove_con0:
    condition0.remove(items)
for items in remove_con1:
    condition1.remove(items)
for items in remove_con2:
    condition2.remove(items)
for items in remove_con3:
    condition3.remove(items)
for items in remove_con4:
    condition4.remove(items)
for items in remove_con5:
    condition5.remove(items)
for items in remove_con6:
    condition6.remove(items)




plt.plot(age,condition2)
plt.show()
