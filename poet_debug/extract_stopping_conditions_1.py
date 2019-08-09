#Extract lines which has stoppping conditions
import os
import numpy


fname='output_poet_debug.txt'
#fname='test_file.txt'
fpath='/mnt/md0/ruskin/QstarFromTidalSynchronization/poet_debug'
out_file=os.path.join(fpath,fname)

fo = open('stopping_condition_output_1.txt','w')
start_name = 'Stored'
stop_name = 'Attempting'
write_line=False
with open(out_file,'r') as f:
    for i,line in enumerate(f):
        if i>20535:
            x=line.split()
            if x:
                if x[0]==start_name:
                    write_line=True
                if x[0]==stop_name:
                    fo.write(line)
                    write_line=False
                if write_line==True:fo.write(line)
