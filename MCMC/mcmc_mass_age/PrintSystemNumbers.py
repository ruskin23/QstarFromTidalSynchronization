s=[]

with open('catalog_KIC.txt','r') as f:
    next(f)
    for lines in f :
        x=lines.split()
        s.append(x[0])

print(s)
