import corner
import numpy as np
import matplotlib.pyplot as plt
import argparse


parser=argparse.ArgumentParser()
parser.add_argument('-l',action='store',dest='system',
                    help='specify the system fom catalog')
args=parser.parse_args()

system=args.system

data=[]

#with open('MassAgeFehSamples_'+system+'.txt','r') as f:
with open('mass_age_teff_sample_8.txt','r') as f:
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
                     labels=[r"$mass$",r"age",r"$feh$",r"$teff$",r"logg"],
                     color='b',
                     qunatiles=[0.16,0.5],
                     show_titles=True)



plt.figure(1)
plt.show()

