import corner
import pickle
import numpy
import matplotlib.pyplot as plt
import sys



# system=sys.argv[1]
s=['1', '8', '12', '13', '17', '20', '25', '28', '32', '36', '39', '43', '44', '47', '48', '50', '54', '56', '57', '67', '70', '73', '76', '79', '80', '81', '83', '84', '85', '86', '88', '92', '93', '94', '95', '96', '106', '109', '120', '123', '126', '137']
s=['50']

with open('complete_chains.pickle','rb') as f:
    D=pickle.load(f)
for system in s:
    data=[]

    print(f'System {system}')
    for system_name,parameters in D.items():
        if system_name==system:
            for param,value in parameters.items():
                if param in ['Porb','logQ','Spin']:
                    data.append(value)

    d=numpy.vstack(data)
    d=d.T

    figure=corner.corner(d,
                        # labels=[r"$Porb$",r"eccentricity",r"$Wdisk$",r"$logQ$",r"$mass$",r"$age$",r"$feh$",r"$Pspin$"],
                        labels=[r"$P_{orb}$",r"$\log_{10}{Q}$",r"$P_{star}$"],
                        color='r',
                        show_titles=True)

    plt.figure(1)
    # plt.title('KIC 6579806')
    # plt.show()
    # plt.savefig('cornor_plots/CornerPlot_'+system+'.png')
    plt.savefig('non_sync.png')
    plt.close()

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
