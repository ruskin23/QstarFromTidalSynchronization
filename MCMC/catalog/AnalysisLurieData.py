with open('Lurie_2017.tsv','r') as f:
    for lines in f:
        x=lines.split()
        if len(x)==5:
            KIC=x[0]
            Porb=float(x[1])
            Pacf=float(x[2])
            P1min=float(x[3])
            P1max=float(x[4])
            r=P1min/Porb
            if abs(r-2)<0.01:
                print(r)
                print(KIC)
