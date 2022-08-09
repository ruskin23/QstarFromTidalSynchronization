from ftplib import MAXLINE
import pandas
from pathlib import Path as path


import numpy
d=[]
with open('/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/version2_emcee/catalog/filtering/nominal_value_catalog_Iconv_cutoff.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        d.append(x[1])

print(len(d))
# print(*d, sep='\t')
for i in range(11):
    print(*d[8*i:8*(i+1)],sep=' ')

