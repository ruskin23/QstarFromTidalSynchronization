import matplotlib.pyplot as plt
import numpy



spin_grid = numpy.empty((73,3))
with open('nominal_value_catalog_Iconv_cutoff.txt', 'r') as f:
    next(f)
    for i, lines in enumerate(f):
        x = lines.split()
        spin_grid[i] = numpy.array([k for k in [x[2],x[3],x[4]]])

spin_grid = numpy.transpose(spin_grid)

spin_error_upper = spin_grid[1] + spin_grid[2]
spin_error_lower = spin_grid[1] - spin_grid[2]

plt.scatter(spin_grid[0], numpy.log10(spin_grid[1]/spin_grid[0]), color = 'r')
plt.plot(spin_grid[0], numpy.zeros(73), color = 'k', linestyle = '--')
plt.xscale('log')
plt.savefig('lurie_2017.png')
        

