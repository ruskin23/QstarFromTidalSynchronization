import re

with open('data_file.txt','w') as fd:
    fd.write('KIC'+'\n')

with open('sol_file.csv','r') as f:
    for line in f:
        x=line.split()
        y=list(x[0])
        z=[]
        for i in range(16,23):
            z.append(y[i])
        KIC = ''.join(z)
        q=x[1]
        with open ('KIC_catalog.txt','r') as f1:
            for line1 in f1:
                if re.search(KIC,line1):
                    data=line1.split()
                    data.append(q)
                    with open('data_file.txt','a') as f2:
                        for d in data:
                            f2.write(d+'\t')
                        f2.write('\n')
                    break

