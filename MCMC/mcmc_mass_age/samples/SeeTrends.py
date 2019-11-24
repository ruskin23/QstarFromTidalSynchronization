import matplotlib.pyplot as plt


age=[]
PrimaryMass=[]
SecondaryMass=[]
feh=[]
eccentricity=[]
Porb=[]
system=[]
with open('SpinlogQCatalog_el0.4.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        system.append(x[0])
        massratio=float(x[14])
        primaryMass=float(x[15])
        if massratio>1:print(x[0])
        secondaryMass=massratio*primaryMass
        PrimaryMass.append(primaryMass)
        SecondaryMass.append(secondaryMass)
        eccentricity.append(float(x[8]))
        Porb.append(float(x[6]))
        age.append(float(x[16]))

for i in range(len(age)):
    if i==0:
        plt.plot(age[i],PrimaryMass[i],color='r',marker='x',label='PrimaryMass')
        plt.plot(age[i],SecondaryMass[i],color='g',marker='x',label='SecondaryMass')
    else:
        plt.plot(age[i],PrimaryMass[i],color='r',marker='x')
        plt.plot(age[i],SecondaryMass[i],color='g',marker='x')
    plt.plot([age[i],age[i]],[PrimaryMass[i],SecondaryMass[i]],color='k')

plt.xlabel('age')
plt.ylabel('Mass')
plt.legend()
plt.savefig('MassVsAge.eps')
