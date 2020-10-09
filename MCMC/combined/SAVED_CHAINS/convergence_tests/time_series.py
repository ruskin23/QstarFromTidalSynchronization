import matplotlib.pyplot as plt


system='76'
clusters=['ganymede','stampede']

parameters=['iteration','porb','eccentricity','wdisk','logQ','mass','age','feh']

for c in clusters:
    for i in range(5):
        n=str(i+1)
        it=[]
        values=[]
        with open('../'+c+'/MCMC_'+system+'/accepted_parameters_'+n+'.txt') as f:
            next(f)
            for lines in f:
                x=lines.split()
                it.append(int(x[0]))
                values.append(float(x[parameters.index('logQ')]))
            plt.scatter(it,values)

plt.show()
