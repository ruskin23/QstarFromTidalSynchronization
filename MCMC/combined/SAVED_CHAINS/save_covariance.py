import shelve
import sys
import pickle
from covariance_matrix import Covariance

#system=sys.argv[1]
covariance=dict()
s=['39', '137', '54', '80', '126', '76', '1', '32', '67' ,'81', '95', '96', '50', '85', '56','73', '86', '92', '20']
for system in s:

    C=Covariance(system).Calculate('Covariance')

    covariance[system]=C

pickle.dump(covariance,open('covariance.pickle','wb'))


#s=shelve.open('covariance.db')

#C=Covariance(system).Calculate('Covariance')
#s[system]=C
#s.close()

