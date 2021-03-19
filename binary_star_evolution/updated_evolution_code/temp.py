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


# import sys
# import os
# from pathlib import Path
# from directories import directories

# home_dir=str(Path.home())
# path=directories(home_dir)
# sys.path.append(path.poet_path+'/PythonPackage')
# sys.path.append(path.poet_path+'/scripts')

# print(f'{path.results_directory}/spin.txt')

import numpy
s=['1', '8', '12', '13', '17', '20', '25', '28', '32', '36', '39', '43', '44', '47', '48', '50', '54', '56', '67', '70', '73', '76', '79', '80', '81', '83', '84', '85', '86', '88', '92', '93', '94', '95', '96', '106', '109', '120', '123', '126', '137']

with open('SpinlogQCatalog_el0.4.txt','r') as f:
    next (f)
    for lines in f:
        x=lines.split()
        if x[0] in s:
            teff_lower=float(x[2])-float(x[3])
            teff_higher=float(x[2])+float(x[3])
            for t in [teff_lower,teff_higher]:
                if numpy.logical_and(t>3500,t<5000):
                    print(f'System = {x[0]} teff_lower  = {teff_lower} teff_higher = {teff_higher}')
                    break