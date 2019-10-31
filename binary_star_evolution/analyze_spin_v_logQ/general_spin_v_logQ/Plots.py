from matplotlib import pyplot as plt
spin=[]
porb=[]
e=[]
r=[]
with open('spin_vs_logQ_systems_0.2.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        spin.append(float(x[12]))
        porb.append(float(x[6]))
        r.append(float(x[12])/float(x[6]))
        e.append(float(x[8]))


plt.scatter(spin,r)
plt.savefig('Ratio.png')
