import corner
import pickle
import numpy
import matplotlib.pyplot as plt
plt.style.use('figureParams.mplstyle')
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

Porb=data[0]
Pspin=data[7]
logQ=data[3]
data1=numpy.array([Porb,logQ,Pspin])
d1=numpy.vstack(data1)
d1=d1.T
print(d1.shape)

figure=corner.corner(d1,
                    labels=[r"$Porb$",r"$\log_{10}{Q^{'}_{*}}$",r"$Pspin$"],
                    color='r',
                    show_titles=True)


# d=numpy.vstack(data)
# d=d.T
# print(d.shape)


# figure=corner.corner(d,
#                     labels=[r"$Porb$",r"eccentricity",r"$Wdisk$",r"$logQ$",r"$mass$",r"$age$",r"$feh$",r"$Pspin$"],
#                     show_titles=True)

plt.savefig('non_sync.pdf')

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
