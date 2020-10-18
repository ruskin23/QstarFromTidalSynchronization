from covariance_kartof import Covariance
import matplotlib.pyplot as plt
from numpy.linalg import inv
import numpy
c=Covariance('39').Calculate('Covariance')

c=numpy.log10(abs(inv(c)))
plt.imshow(c)
plt.colorbar()


plt.show()
