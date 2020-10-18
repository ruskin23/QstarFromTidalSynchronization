from covariance_matrix import Covariance
import numpy


system = '137'

C=Covariance(system)
print(C.Calculate('Covariance'))
print(numpy.linalg.inv(C.Calculate('Covariance')))
