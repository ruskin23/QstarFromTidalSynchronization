import matplotlib.pyplot as plt
import numpy
import sys
import itertools

with open('Ages.txt','w') as f:
    f.write('System'+'\t'+
            '1'+'\t'+
            '2'+'\t'+
            '3'+'\t'+
            '4'+'\t'+
            '5'+'\t'+
            '10'+'\t'+
            '20'+'\t'+
            '30'+'\t'+
            '40'+'\t'+
            '50'+'\n')

def _cummulative_distribution(x):
    a=[]
    v=0.0
    values=[]

    for tuples in x:
        a=numpy.append(a,tuples[0])
        v=v+tuples[1]
        values=numpy.append(values,v)

    return list(zip(a,values/max(values)))

systems=['2', '5', '14', '18', '23', '27', '33', '35', '37', '49', '59', '62', '63', '65', '68', '69', '77', '78', '82', '99', '100', '111', '117', '119', '125', '128', '129', '132', '133', '134', '138', '142']
#systems=['23']
for system in systems:
    filename='updated_samples/MassAgeFehSamples_'+system+'.txt'

    ages=[]
    multiplicity=[]
    with open(filename,'r') as f:
        next(f)
        for lines in f:
            x=lines.split()
            ages=numpy.append(ages,float(x[1]))
            multiplicity=numpy.append(multiplicity,float(x[3]))


    ages_tuples=zip(ages,multiplicity)
    ages_tuples=sorted(ages_tuples, key=lambda tup: tup[0])

    ages_cummulative=_cummulative_distribution(ages_tuples)
    percentile=[0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5]
    with open('Ages.txt','a') as f:
        f.write(system+'\t')
        a=[]
        for p in percentile:
            for age in ages_cummulative:
                if abs(p-age[1])<0.01:
                    a=numpy.append(a,str(age[0]))
                    break
        f.write('\t'.join(a))
        f.write('\n')
