import shelve


git_dir='/QstarFromTidalSynchronization/MCMC/combined'
work_dir='/work/06850/rpatel23/stampede2'

current_directory=work_dir+git_dir

s=shelve.open('covariance.db')
covariance_matrix=s['81']
s.close()

print(covariance_matrix)
