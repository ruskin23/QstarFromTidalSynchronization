import re
import sys
import os
import argparse

f2=open('cross_2.txt','w')
with open('cross_1.txt') as f:
    for line in f:
        x=line.split()
        if len(x)==5:
            f2.write(line)
f2.close()
