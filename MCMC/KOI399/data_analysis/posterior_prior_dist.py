import scipy
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import sys

instance = sys.argv[1]

distribution = dict(
                        age=dict(value=4.6, sigma=3.0),
                        teff_primary=dict(value=5922.0, sigma=200.0),
                        feh=dict(value=-0.06, sigma=0.11),
                        Wdisk=dict(value=2*scipy.pi / 1.4, sigma=0.1),
                        logQ = dict(value=0.0,sigma=0.0)
                    )

name = 'accepted_parameters_'+instance+'.txt'

all_array = []
for index,key in enumerate(distribution):
    array=[]
    index = index+1
    with open(name, 'r') as f:
        next(f)
        for line in f:
            row = line.split()
            array.append(float(row[index]))
    all_array.append(array)
    distribution[key]['value']= np.mean(array)
    distribution[key]['sigma'] = np.std(array)

prior = dict(
                        age=dict(value=4.6, sigma=3.0),
                        teff_primary=dict(value=5922.0, sigma=200.0),
                        feh=dict(value=-0.06, sigma=0.11),
                        Wdisk=dict(value=2*scipy.pi / 1.4, sigma=0.1)
                    )


teff_range = np.linspace(5322,6522,1000)
t_post = norm.pdf(teff_range,loc = distribution['teff_primary']['value'],scale
                   = distribution['teff_primary']['sigma'])
t_prior = norm.pdf(teff_range,loc = prior['teff_primary']['value'],scale
                   = prior['teff_primary']['sigma'])




plt.plot(teff_range,t_post,'-g')
plt.plot(teff_range,t_prior,'-r')
plt.show()


#plt.hist2d(all_array[2],all_array[4],bins=10)
#cbar = plt.colorbar()
#cbar.ax.set_ylabel('Counts')
#plt.show()
