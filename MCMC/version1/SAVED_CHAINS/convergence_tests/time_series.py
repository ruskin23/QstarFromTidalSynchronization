import matplotlib.pyplot as plt
from ks_tests import _fill_parameters


system='84'
clusters=['ganymede','stampede']

parameters=['iteration','porb','eccentricity','wdisk','logQ','mass','age','feh']

for c in clusters:
    for i in range(5):
        n=str(i+1)
        it=[]
        values=[]
        filename='../'+c+'/MCMC_'+system+'/accepted_parameters_'+n+'.txt'
        filled_file=_fill_parameters(filename)
        with open(filled_file,'r') as f:
            #next(f)
            for lines in f:
                x=lines.split()
                it.append(int(x[0]))
                values.append(float(x[parameters.index('logQ')]))
            plt.plot(it,values)

            plt.show()
