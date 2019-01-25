import numpy
import os
import re
import csv

filename = 'test_3.txt'


result = 'A'
x = numpy.linspace(1,100,5)

os.remove('test_3.txt')
for i in range(1,2):
    with open('test_3.txt','a') as file:
        for n in x:
            file.write(repr(n*i) + '\t')
        file.write('\n')
    file.close()



with open('test_3.txt', 'r') as f:
        
    reader = csv.reader(f, dialect='excel-tab')
    for row in reader:
        array = row

f.close()



import argparse

parser = argparse.ArgumentParser()

parser.add_argument("a", nargs='?', default="check_string_for_empty")

args = parser.parse_args()

if args.a == 'check_string_for_empty':
    print('default') 
elif args.a == 'start':
    print('begin') 
else:
    print (args.a)



