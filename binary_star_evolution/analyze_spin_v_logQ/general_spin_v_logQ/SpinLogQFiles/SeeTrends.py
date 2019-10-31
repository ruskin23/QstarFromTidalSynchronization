import numpy
import matplotlib.pyplot as plt

system_array=[]
PspinValue=[]
with open('spin_vs_logQ_systems_0.2.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        system_array.append(x[0])
        PspinValue.append(float(x[12]))


no_solution=[]
one_solution=[]
multiple_solution=[]
for spin,system in zip(PspinValue,system_array):
    SpinFilename='SpinLogQ_'+system+'.txt'
    with open(SpinFilename,'r') as f1:
        next(f1)
        p=[]
        for i,lines in enumerate(f1):
            x=lines.split()
            p.append(float(x[1])-spin)
        p=numpy.array(p)
        zero_crossing=numpy.where(numpy.diff(numpy.sign(p)))[0]
        if zero_crossing.size==0:
            no_solution.append(system)

        elif zero_crossing.size==1:
            one_solution.append(system)

        elif zero_crossing.size>1:
            multiple_solution.append(system)

print(no_solution)
print(one_solution)
print(multiple_solution)
#spin,prob ratio

def get_ratios(solutions,flag):
    ratios=[]
    eccentrcity=[]
    for system in solutions:
        with open('spin_vs_logQ_systems_0.2.txt','r') as f:
            for lines in f:
                x=lines.split()
                if x[0]==system:
                    r=float(x[6])/float(x[12])
                    ratios.append(r)
                    eccentrcity.append(float(x[8]))
    if flag=='r':return ratios
    if flag=='e':return eccentrcity

no_solution_ratio=get_ratios(no_solution,'r')
one_solution_ratio=get_ratios(one_solution,'r')
multiple_solution_ratio=get_ratios(multiple_solution,'r')
print(len(no_solution))
print(len(one_solution))
print(len(multiple_solution))

#eccentrcity

no_solution_eccentricty=get_ratios(no_solution,'e')
one_solution_eccentricty=get_ratios(one_solution,'e')
multiple_solution_eccentricty=get_ratios(multiple_solution,'e')


no_x=len(no_solution)*[1]
one_x=len(one_solution)*[1]
multi_x=len(multiple_solution)*[1]

plt.scatter(no_solution_ratio,no_x,color='r')
plt.scatter(one_solution_ratio,one_x,color='g')
plt.scatter(multiple_solution_ratio,multi_x,color='b')
plt.show()

plt.scatter(no_solution_eccentricty,no_x,color='r')
plt.scatter(one_solution_eccentricty,one_x,one_x,color='g')
plt.scatter(multiple_solution_eccentricty,multi_x,multi_x,color='b')
plt.show()

