import corner
import numpy as np
import matplotlib.pyplot as plt

data=[]
with open('AccetedParameters.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        i=1
        data.append(float(x[1]))	
        data.append(float(x[2]))
        data.append(float(x[3]))
        data.append(float(x[4]))
        data.append(float(x[5]))
        data.append(float(x[6]))
        data.append(float(x[7]))
        data.append(float(x[13]))

data=np.array(data)
d=data.reshape([len(data)//8,8])

figure=corner.corner(d,
                     labels=[r"$mass$",r"age",r"$feh$",r"$Porb$",r"$e$",r"$Wdisk$",r"$logQ$",r"$Pspin$"],
                     color='k',
                     show_titles=True)

plt.figure(1)
plt.show()
