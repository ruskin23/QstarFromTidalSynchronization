import matplotlib.pyplot as plt

p=[]
e=[]

with open('CommonKICandWindeCat.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        p.append(float(x[2]))
        e.append(float(x[4]))

plt.scatter(p,e)
#plt.yscale('log')
plt.show()
