import numpy

p=[]
with open('extracted.txt','r') as f:
    for lines in f:
        x=lines.split()

        p=numpy.append(p,float(x[0]))


print(numpy.mean(p))


