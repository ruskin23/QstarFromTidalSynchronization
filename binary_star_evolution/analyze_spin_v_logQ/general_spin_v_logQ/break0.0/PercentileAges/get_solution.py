import numpy
import sys
sys.path.append('../../../../../MCMC/mcmc_mass_age/samples/')
from PercentileClass import PercentileAge
import matplotlib.pyplot as plt


system=sys.argv[1]
percentiles=sys.argv[2]
percentiles=percentiles.split(',')

Age=PercentileAge(system)

solution_file='Solutions_'+system+'.txt'
with open(solution_file,'w') as f:
    f.write('percentile'+'\t'+
            'age'+'\t'+
            'logQ'+'\t'+
            'Pspin'+'\t'+
            'error'+'\n')

def observed_values():

    with open('../../SpinlogQCatalog_el0.4.txt','r') as f:
        for i,lines in enumerate(f):
            if i==0:
                x=lines.split()
                spin_index=x.index('Pspin')
                orbital_index=x.index('Porb')
            else:
                x=lines.split()
                if x[0]==system:
                    Pspin=float(x[spin_index])
                    PspinError=float(x[spin_index+1])
                    Porb=float(x[orbital_index])
                    break

    return Pspin



def get_solution(percentile,age):

    spinvlogQfile='System_'+system+'/SpinLogQ_'+str(percentile)+'.txt'
    q=[]
    p=[]
    with open(spinvlogQfile,'r') as f:
        next(f)
        for lines in f:
            x=lines.split()
            if numpy.logical_and(abs(float(x[6]))<1e-3,
                                 abs(float(x[7]))<1e-3):
                q=numpy.append(q,float(x[0]))
                p=numpy.append(p,float(x[1]))

    data=zip(q,p)
    data=sorted(data,key=lambda tup: tup[0])
    q_values=[]
    p_values=[]
    for t in data:
        q_values=numpy.append(q_values,t[0])
        p_values=numpy.append(p_values,t[1])
    qmin=min(q_values)
    qmax=max(q_values)
    print('For age {} q0={} qn={}'.format(percentile,qmin,qmax))

    q_array=numpy.linspace(qmin,qmax,10000)
    p_interpolated=numpy.interp(q_array,q_values,p_values)

    p_observed=observed_values()
    p_diff=numpy.array(p_interpolated)-p_observed

    zero_crossing=numpy.where(numpy.diff(numpy.sign(p_diff)))[0]
    sol_p=p_interpolated[zero_crossing][0]
    sol_q=q_array[zero_crossing][0]
    error_sol_p=abs(sol_p-p_observed)

    with open(solution_file,'a') as f:
        f.write(repr(percentile)+'\t'+
                repr(age)+'\t'+
                repr(sol_q)+'\t'+
                repr(sol_p)+'\t'+
                repr(error_sol_p)+'\n')

for p in percentiles:
    percentile=int(p)
    age=Age(percentile)
    get_solution(percentile,age)

