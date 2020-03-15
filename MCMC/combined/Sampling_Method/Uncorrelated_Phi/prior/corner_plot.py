import corner
import numpy as np
import matplotlib.pyplot as plt

data=[]

logQ=[]
age=[]
with open('AcceptedParameters.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        data.append(float(x[1]))
        data.append(float(x[2]))
        data.append(float(x[3]))
        data.append(float(x[4]))
        data.append(float(x[5]))
        data.append(float(x[6]))
        data.append(float(x[7]))

data=np.array(data)
d=data.reshape([len(data)//7,7])

figure=corner.corner(d,
                     labels=[r"$Porb$",r"eccentricity",r"$Wdisk$",r"$logQ$",r"$mass$",r"$age$",r"$feh$"],
                     color='k',
                     show_titles=True)



plt.figure(1)
plt.show()

