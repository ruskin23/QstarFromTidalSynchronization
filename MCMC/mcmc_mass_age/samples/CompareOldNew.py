import matplotlib.pyplot as plt


OldAge=[]
NewAge=[]
systemNew=[]
systemOld=[]

with open('SpinlogQCatalog_el0.4.txt','r') as fN:
    next(fN)
    for linesNew in fN:
        x=linesNew.split()
        with open('OldCatalog.txt','r') as fO:
            next(fO)
            for linesOld in fO:
                y=linesOld.split()
                if x[0]==y[0]:
                    print('Found System = ',systemNew)
                    systemNew.append(int(x[0]))
                    systemOld.append(int(y[0]))
                    NewAge.append(float(x[17]))
                    OldAge.append(float(y[17]))
                    break


for i in range(len(systemNew)):
    if i==0:
        plt.plot(systemOld[i],OldAge[i],color='r',marker='x',label='OldSystem')
        plt.plot(systemNew[i],NewAge[i],color='g',marker='x',label='NewSystem')
    else:
        plt.plot(systemOld[i],OldAge[i],color='r',marker='x')
        plt.plot(systemNew[i],NewAge[i],color='g',marker='x')
    plt.plot([systemOld[i],systemNew[i]],[OldAge[i],NewAge[i]],color='k')
plt.xlabel('System Number')
plt.ylabel('FeH')
plt.legend()
#plt.show()
plt.savefig('CompareFeH.eps')
