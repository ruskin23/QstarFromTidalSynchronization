import argparse
import sys

key = {'a':['a','b','c'],
       'b':[4,5,6],
       'c':[7,8,9]}

import csv
with open('text.csv', 'w') as f:
    writer = csv.writer(f, delimiter='\t')
#    for k,array in key.items():
    writer.writerows(zip(*key.values()))

