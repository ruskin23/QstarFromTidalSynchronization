#Extra formatting required in stopping condition
import os
import numpy


fname_in='stopping_condition_output_2.txt'
fname_out='stopping_condition_output_3.txt'
fpath='/mnt/md0/ruskin/QstarFromTidalSynchronization/poet_debug'


in_file=os.path.join(fpath,fname_in)
out_file=os.path.join(fpath,fname_out)

fo=open(out_file,'w')


itmes=['Age:|','Age:','Condition[','Derivative[','NOT','Condition','SELECTED']

with open(in_file,'r') as f:
    for line in f:
        d=line.split()
        if d:
            if d[0] in itmes:
                fo.write(line)



fo.close()



