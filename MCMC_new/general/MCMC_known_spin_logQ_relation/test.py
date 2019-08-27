with open('sol_file.csv','r') as f:
    for line in f:
        x=line.split()
        y=list(x[0])
        z=[]
        for i in range(16,23):
            z.append(y[i])
        KIC = ''.join(z)
        print(KIC)
