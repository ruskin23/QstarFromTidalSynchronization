

viable=True
with open('NewCatalog.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        if float(x[2])>30:viable=False
        if float(x[4])>0.45:viable=False
        if float(x[6])>1.2 and float(x[6])<0.4:viable=False
        if float(x[8])>1.2 and float(x[8])<0.4:viable=False
        if viable==True:
            print(x[0])
        else:viable=True

