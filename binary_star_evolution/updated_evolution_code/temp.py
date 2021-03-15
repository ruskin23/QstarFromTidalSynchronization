# import numpy

# with open('params.txt','w') as ff:
#     ff.write('porb_initial\te_initial\tWdisk\tlogQ\tprimary_mass\tage\tfeh\n')
#     with open('temp.out','r') as f:
#         for lines in f:
#             p=lines.split()
#             if float(p[5][:-1]) > 4*numpy.pi/float(p[1][:-1]):cond_satisfied=True
#             else:cond_satisfied=False
#             ff.write(repr(float(p[1][:-1]))+'\t'+
#                      repr(float(p[3][:-1]))+'\t'+
#                      repr(float(p[5][:-1]))+'\t'+
#                      repr(float(p[7][:-1]))+'\t'+
#                      repr(float(p[9][:-1]))+'\t'+
#                      repr(float(p[11][:-1]))+'\t'+
#                      repr(float(p[13][:-1]))+'\t'+
#                      repr(cond_satisfied)+'\n'
#                         )


import sys
import os
from pathlib import Path
from directories import directories

home_dir=str(Path.home())
path=directories(home_dir)
sys.path.append(path.poet_path+'/PythonPackage')
sys.path.append(path.poet_path+'/scripts')

print(f'{path.results_directory}/spin.txt')