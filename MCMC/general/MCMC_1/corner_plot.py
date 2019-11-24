import corner
import numpy as np
import matplotlib.pyplot as plt

data=[]
with open('AccetedParameters.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        i=1
        while i<8:
            data.append(float(x[i]))
            i=i+1

data=np.array(data)
d=data.reshape([len(data)//7,7])

figure=corner.corner(d,
                     labels=[r"$mass$",r"age",r"$feh$",r"$Porb$",r"$e$",r"$Wdisk$",r"$logQ$"],
                     color='k',
                     show_titles=True)

plt.figure(1)
plt.show()
