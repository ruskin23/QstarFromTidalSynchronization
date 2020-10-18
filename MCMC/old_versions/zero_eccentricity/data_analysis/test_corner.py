import corner
import numpy as np
import matplotlib.pyplot as plt

data1 = []
data2 = []

with open('combined_accepted_parameters.txt','r') as f:
    for lines in f:
        data = lines.split()
        data1.append(float(data[1]))
        data1.append(float(data[2]))
        data1.append(float(data[3]))
        data1.append(float(data[4]))
        data1.append(float(data[5]))
        data1.append(float(data[6]))

data1=np.array(data1)
data = data1.reshape([len(data1)//6,6])
#d = np.vstack([data1,data2])
figure = corner.corner(data,
                       labels=[ r"$Teff$", r"$FeH$", r"$Porb", r"$logg", r"$Wdisk$",r"$\log Q$"],
                       color='b',
                       show_titles=True)
plt.figure(1)

figure.savefig('plot.png')
plt.show()
