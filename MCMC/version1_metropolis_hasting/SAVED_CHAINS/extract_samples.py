import pickle


parameters=['Porb','eccentricity','Wdisk','logQ','primary_mass','age','feh','Spin']

with open('complete_chains.pickle','rb') as f:
    D=pickle.load(f)

data=[]

for system_name,parameters in D.items():
    if system_name=='1':
        for param,value in parameters.items():
            data.append(value)

print(data[0])