import shelve
import sys
from covariance_matrix import Covariance

system=sys.argv[1]

s=shelve.open('covariance.db')

C=Covariance(system).Calculate('Covariance')
print(C)
s[system]=C
s.close()

