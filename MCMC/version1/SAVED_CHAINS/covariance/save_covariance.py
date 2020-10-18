import shelve
import sys
import pickle
from covariance_matrix import Covariance

#system=sys.argv[1]

#Create covariance matrix for first time
#covariance=dict()
s=['39', '137', '54', '80', '126', '76', '1', '32', '67' ,'81', '95', '96', '50', '85', '56','73', '86', '92', '20']
#for system in s:
#
#    C=Covariance(system).Calculate('Covariance')
#
#    covariance[system]=C
#
#pickle.dump(covariance,open('covariance.pickle','wb'))


#Update covariance matrix for 76 54 39 92
updated_covariance=dict()
no_update=['137', '80', '126', '1', '32', '67' ,'81', '95', '96', '50', '85', '56','73', '86', '20']
update=['39', '54', '76', '92']
for system in s:
    if system in no_update:
        updated_covariance[system]=pickle.load(open('covariance.pickle','rb'))[system]
    else:
        updated_covariance[system]=Covariance(system).Calculate('Covariance')
pickle.dump(updated_covariance,open('updated_covariance.pickle','wb'))

