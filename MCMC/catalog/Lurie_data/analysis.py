import numpy


n_out=0
n_in=0
with open('filter2.txt','r') as f:
    for lines in f:
        x=lines.split()
        porb=float(x[1])
        pspin=float(x[5])
        worb=2*numpy.pi/porb
        wspin=2*numpy.pi/pspin
        r=worb/wspin
        if r<2:n_in=n_in+1
        else:n_out=n_out+1
        

print(f'in:{n_in} out:{n_out}')
