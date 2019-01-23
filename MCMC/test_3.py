import numpy



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
