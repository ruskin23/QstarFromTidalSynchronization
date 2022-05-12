import sys
import os
import random
import numpy
sys.path.append('/home/ruskin/projects/poet/PythonPackage')
sys.path.append('/home/ruskin/projects/poet/scripts')

from stellar_evolution.library_interface import library

from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as \
    orbital_evolution_library

#Put Porb<60 days filter from raw orbital catalog
fo=open('filter_porb.txt','w')
with open('windemuth_orbital_raw.txt','r') as f:
    for i,lines in enumerate(f):
        x=lines.split()
        if i<2:fo.write(lines)
        else:
            if float(x[1])<60:fo.write(lines)
fo.close()

#Combine filtered Porb<60 days catalog with stellar parameters catalog
fo=open('combined_porb_filtered_catalog.txt','w')
fo.write('Number'+'\t'+
         'KIC'+'\t'+
         'Porb'+'\t'+
         'esinw'+'\t'+
         'esinw+sigma'+'\t'+
         'esinw-sigma'+'\t'+
         'ecosw'+'\t'+
         'ecos+sigma'+'\t'+
         'ecos-sigma'+'\t'+
         'z'+'\t'+
         'z+sigma'+'\t'+
         'z-sigma'+'\t'+
         'age'+'\t'+
         'age+sigma'+'\t'+
         'age-sigma'+'\t'+
         'm1+sigma'+'\t'+
         'm1-sigma'+'\t'+
         'm2+sigma'+'\t'+
         'm2-sigma'+'\n')

match=False
with open('filter_porb.txt','r') as fp:
    next(fp)
    next(fp)
    for i,lines1 in enumerate(fp):
        params=[]
        x=lines1.split()
        params.append(str(i+1))
        params.append(x[0])
        params.append(x[1])
        params.append(x[7])
        params.append(x[8])
        params.append(x[9])
        params.append(x[10])
        params.append(x[11])
        params.append(x[12])
        with open('windemuth_stellar_raw.txt','r') as fs:
            next(fs)
            next(fs)
            for lines2 in fs:
                y=lines2.split()
                if y[0]==x[0]:
                    params.append(y[1])
                    params.append(y[2])
                    params.append(y[3])
                    params.append(y[4])
                    params.append(y[5])
                    params.append(y[6])
                    params.append(y[7])
                    params.append(y[8])
                    params.append(y[9])
                    params.append(y[10])
                    params.append(y[11])
                    params.append(y[12])
                    line='\t'.join(params)
                    fo.write(line)
                    fo.write('\n')


#catalog of nominal values
#condition: Porb<60
fo=open('normandy_catalog.txt','w')
fo.write('Number'+'\t'+
         'KIC'+'\t'+
         'Porb'+'\t'+
         'eccentricity'+'\t'+
         'z'+'\t'+
         'age'+'\t'+
         'm1'+'\t'+
         'm2'+'\n')

with open('combined_porb_filtered_catalog.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        esinw=float(x[3])
        ecosw=float(x[6])
        eccentricity=numpy.sqrt(esinw**2 + ecosw**2)
        if eccentricity<0.50 and float(x[12])>8:
            fo.write(x[0]+'\t'+x[1]+'\t'+x[2]+'\t'+repr(eccentricity)+'\t'+x[9]+'\t'+x[12]+'\t'+x[15]+'\t'+x[18]+'\n')


#catalog with following conditions:
#Porb<60
#m1 and m2 in POET range
#feh in poet range
#condition: Porb<60
fo=open('normandy_poet_catalog.txt','w')
fo.write('Number'+'\t'+
         'KIC'+'\t'+
         'Porb'+'\t'+
         'eccentricity'+'\t'+
         'z'+'\t'+
         'age'+'\t'+
         'm1'+'\t'+
         'm2'+'\n')
KIC=[]
with open('normandy_catalog.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        m1=float(x[6])
        m2=float(x[7])
        z=float(x[4])
        feh=library.feh_from_z(z)
        if  numpy.logical_and(
                            numpy.logical_and(feh>-1.014,feh<0.537)
                            ,
                            numpy.logical_and(numpy.logical_and(m1>0.4,m1<1.2),
                                              numpy.logical_and(m2>0.4,m2<1.2))
        )                   :
            KIC.append(x[1])
            fo.write(lines)


print(len(KIC))
f_final=open('spin_catalog.txt','w')
n=0
with open('Lurie_binaries_with_p1.txt','r') as f:
    for lines in f:
        x=lines.split()
        if x[0] in KIC:
            n=n+1
            spin=float(x[5])
            spin_error=abs(float(x[6])-float(x[5]))
print(n)
