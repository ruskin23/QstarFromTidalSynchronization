import sys
import numpy
import scipy
from scipy import optimize

def function(M2,M1,f_M,i):
    return (((M2*numpy.sin(i))**1.5)/numpy.sqrt(f_M)) - M2 - M1


def mass_limits(M1,f_M):

    M=numpy.linspace(0.4,1.6,10)
    previous_guess=0.4
    previous_value=function(previous_guess,M1,f_M)

    for m in M:

        currect_value=function(m,M1,f_M)
        if previous_value*currect_value<0:
            return previous_guess,m
        else:
            previous_value=currect_value
            previous_guess=m
            continue


def secondary(M1,f_M,i):


    #M2_i,M2_f=mass_limits(M1,f_M)
    try:
        sol=scipy.optimize.root_scalar(lambda M2: function(M2,M1,f_M,i),bracket=[0.4,1.2],method='brentq')

        return sol.root
    except ValueError:
        return scipy.nan




#cos_i=numpy.linspace(-1,1,10)
inclination=numpy.linspace(0,numpy.pi,100)
system=sys.argv[1]
with open(system+'_mass_table.txt','w') as f:
    f.write('index'+'\t'+
            'cosi'+'\t'+
            'i'+'\t'+
            'mass_function'+'\t'+
            'primary_mass'+'\t'+
            'secondary_mass'+'\n')

k=0
for i in inclination:
    with open('rv_fit.txt','r') as f:
        next(f)
        for lines in f:
            x=lines.split()
            if x[0]==system:
                K=float(x[1])
                e=float(x[2])
                with open('Table1.txt','r') as f1:
                    next(f1)
                    for lines in f1:
                        y=lines.split()
                        if y[0]==x[0]:
                            Porb=float(y[4])
                            f_M=Porb*K*K*K*((1-e*e)**1.5)*1.03e-7
                            M1=float(y[9])
                            M2=secondary(M1,f_M,i)
                            if numpy.isnan(M2): continue
                            else:
                                k=k+1
                                with open(system+'_mass_table.txt','a') as fw:
                                    fw.write(repr(k)+'\t'+
                                            repr(numpy.cos(i))+'\t'+
                                            repr(i)+'\t'+
                                            repr(f_M)+'\t'+
                                            repr(M1)+'\t'+
                                            repr(M2)+'\n')
