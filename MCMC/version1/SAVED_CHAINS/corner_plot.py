import corner
import pickle
import numpy
import matplotlib.pyplot as plt
import sys



system=sys.argv[1]

with open('complete_chains.pickle','rb') as f:
    D=pickle.load(f)
data=[]

print(f'System {system}')
for system_name,parameters in D.items():
    if system_name==system:
        for param,value in parameters.items():
            data.append(value)

d=numpy.vstack(data)
d=d.T

figure=corner.corner(d,
                    labels=[r"$Porb$",r"eccentricity",r"$Wdisk$",r"$logQ$",r"$mass$",r"$age$",r"$feh$",r"$Pspin$"],
                    show_titles=True)

plt.show()

# with open('AcceptedParameters.txt','r') as f:
#     next(f)
#     for lines in f:
#         x=lines.split()
#         data.append(float(x[1]))
#         data.append(float(x[2]))
#         data.append(float(x[3]))
#         data.append(float(x[4]))
#         data.append(float(x[5]))
#         data.append(float(x[6]))
#         data.append(float(x[7]))
#         data.append(float(x[15]))

# data=np.array(data)
# d=data.reshape([len(data)//8,8])
