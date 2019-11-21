import os
import numpy

index=[]
with open('catalog_KIC.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        index.append(x[0])


print(index)


for i in range(110,113):
    sample_file='mass_age_teff_sample_'+str(i)+'.txt'
    real_index=index[i-1]
    new_file_name='MassAgeFehSamples_'+real_index+'.txt'
    os.rename(sample_file,new_file_name)

