import sys
import pickle
from covariance_matrix import Covariance
import numpy

print(pickle.load(open('covariance.pickle','rb'))['76'])
print(pickle.load(open('updated_covariance.pickle','rb'))['76'])
