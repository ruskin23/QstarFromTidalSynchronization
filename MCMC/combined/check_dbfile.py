import shelve




s=shelve.open('covariance.db')
covariance_matrix=s['81']
s.close()

print(covariance_matrix)
