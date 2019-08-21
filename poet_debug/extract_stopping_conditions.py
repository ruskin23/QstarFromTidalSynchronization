#Extract lines which has stoppping conditions
import os
import numpy


fname='output_poet_debug.txt'
#fname='test_file.txt'
fpath='/mnt/md0/ruskin/QstarFromTidalSynchronization/poet_debug'
out_file=os.path.join(fpath,fname)

fo = open('output.txt','w')
start_name = 'Stored'
stop_name = 'Stepped'
write_line=False
with open(out_file,'r') as f:
    for line in f:
        x=line.split()
        if x:
            if x[0]==start_name:
                write_line=True
            if x[0]==stop_name:
                fo.write(line)
                write_line=False
            if write_line==True:fo.write(line)
