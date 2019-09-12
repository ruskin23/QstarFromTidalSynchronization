import csv
import numpy as np
import scipy
from scipy.stats import norm
import argparse
import sys
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
parser.add_argument('fname',help='enter name of the file')
parser.add_argument('pname',help='enter the name of parameter')
args = parser.parse_args()



parameter = [ 'teff', 'feh', 'Porb', 'logg' ,'wdisk', 'logQ', 'pspin']
try: p = parameter.index(args.pname)
except: print('give correct parameter name')

value =[]

with open(args.fname,'r') as f:
    reader = csv.reader(f, dialect='excel-tab')
    for line in reader:
        value.append(float(line[p+1]))

value.sort()

fig, ax = plt.subplots(figsize=(8, 4))

#histogram of data
n_bins = 100
n, bins, patches = ax.hist(value, bins=n_bins, density=True, cumulative=True, histtype='step')

#cummulative distribution from array
mu = np.mean(value)
sigma = np.std(value)

textstr = '\n'.join( ('$\mu = %f$'%mu,'$\sigma = %f$'%sigma ))
ax.text(0.05,0.95,textstr,transform=ax.transAxes, fontsize=14,
        verticalalignment='top')

y = ((1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-0.5 * (1 / sigma * (bins - mu))**2))
y = y.cumsum()
y /= y[-1]

ax.plot(bins, y, 'k--', linewidth=1.5)



ax.grid(True)
ax.legend(loc='right')
ax.set_title('Cumulative step histograms')
ax.set_xlabel(args.pname)
ax.set_ylabel('Likelihood of occurrence')

plt.show()




