import matplotlib.pyplot as plt

param=2

sample_file_name='AccetedParameters.txt'
multiplicity=0
N=-1
DATA=[]

with open(sample_file_name,'r') as f:
    next(f)
    for i,lines in enumerate(f):
        x=lines.split()
        current_state=[float(x[k]) for k in range(1,5)]
        if i==0:previous_state=current_state
        if previous_state==current_state:
            multiplicity=multiplicity+1
        else:
            previous_state.append(multiplicity)
            DATA.append(previous_state)
            previous_state=current_state
            multiplicity=1
        N=N+1

sorted_DATA=sorted(DATA,key = lambda d: d[param])

value=[]
m=[]
n=0
for d in sorted_DATA:
    value.append(d[param])
    n=n+d[4]
    m.append(n/N)
    print(d)

plt.plot(value,m)
plt.show()
