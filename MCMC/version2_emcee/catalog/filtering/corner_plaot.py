from turtle import shape
import numpy
import corner
import matplotlib.pyplot as plt


A=numpy.load('/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/version2_emcee/catalog/samples/chains/7838639.npz')
all_chains=numpy.transpose(A['thinned_chain'])

esinw_samples=all_chains[8]
ecosw_samples=all_chains[9]


d=numpy.array([esinw_samples,ecosw_samples])
d=numpy.transpose(d)
print(numpy.mean(ecosw_samples))
print(numpy.mean(esinw_samples))
figure=corner.corner(d)
plt.savefig('test.png')