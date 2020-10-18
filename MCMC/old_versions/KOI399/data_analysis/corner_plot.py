import corner
import numpy as np
import matplotlib.pyplot as plt
import argparse



data=[]

with open('combined_accepted_parameters.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        try:
            data.append(float(x[1]))
            data.append(float(x[2]))
            data.append(float(x[3]))
            data.append(float(x[4]))
            data.append(float(x[5]))
        except:continue


data=np.array(data)
d=data.reshape([len(data)//5,5])

figure=corner.corner(d,
                     labels=[r"$age$",r"$teff$",r"feh",r"$Wdisk$",r"$logQ$"],
                     color='b',
                     #qunatiles=[0.16,0.5],
                     show_titles=True)



plt.figure(1)
plt.show()

