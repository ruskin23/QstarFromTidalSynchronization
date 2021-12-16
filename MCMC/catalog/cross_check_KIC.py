import re
import sys
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('fname1',help='enter the name of first file')
parser.add_argument('fname2',help='enter the name of second file')
args = parser.parse_args()

fname1=args.fname1
fname2=args.fname2

f3=open('cross_1_test.txt','w')
with open(fname1,'r') as f1:
    for line in f1:
        name = line.split()
        kic=name[0]
        with open(fname2,'r') as f2:
            for line2 in f2:
                if re.search(kic,line2):
                    #print(line2)
                    f3.write(line2)
                    break
                else:continue
            #    for line2 in f2:
        #        if re.search(kic,line):
        #            f3.write(line2)

f3.close()

if os.stat('cross_1_test.txt').st_size != 0:
    with open('cross_1_test.txt','r') as f:
        for line in f:
            print(line)
