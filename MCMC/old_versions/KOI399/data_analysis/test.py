import csv
import numpy as np
import scipy
from scipy.stats import norm
import argparse
import sys
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

np.seterr(divide='ignore', invalid='ignore')

def func(bins, mu, sigma):
    y = ((1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-0.5 * (1 / sigma * (bins - mu))**2))
    y = y.cumsum()
    y /= y[-1]
    return y


parser = argparse.ArgumentParser()
parser.add_argument('fname',help='enter name of the file')
parser.add_argument('pname',help='enter the name of parameter')
args = parser.parse_args()



parameter = ['age', 'teff', 'feh', 'wdisk', 'logQ', 'pspin']
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

bins = np.delete(bins,100)
xdata=bins
ydata = n

popt, pcov = curve_fit(func, xdata, ydata)

mu = popt[0]
sigma = popt[1]

textstr = '\n'.join( ('$\mu = %f$'%mu,'$\sigma = %f$'%sigma ))
ax.text(0.05,0.95,textstr,transform=ax.transAxes, fontsize=14,
        verticalalignment='top')

ax.plot(bins, func(bins,mu,sigma), 'k--', linewidth=1.5)

ax.grid(True)
ax.set_xlabel(args.pname)
plt.show()
