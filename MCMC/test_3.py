import numpy
import os
import re

filename = 'test_3.txt'


result = 'A'
x = numpy.linspace(10,20,10)

for i in range(x.size):

    with open(filename,'a') as file : 
        
        file.write(repr(i) + "\t" )
        for j in x:
            file.write('%s\t' % j)
        file.write('\t'+result + '\n')

    file.close()



f = os.getcwd()
file_list = os.listdir(f)
name = ["t1.txt" , "t2.txt"]
for n in name:
    if  n in file_list:
        os.remove(n)



import scipy


def return_nan():
    return scipy.nan

check_array = [2 , 3,5 , return_nan(), 12, 15]

for x in check_array:
    if numpy.isnan(x):    
        print("\nnan found")
        continue
    print("\n", x)
