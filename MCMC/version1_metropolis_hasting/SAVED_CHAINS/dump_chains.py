#this code creates a pickle file with chains with filled values for all parameters
#excluding the iterations for which there were consistent creahsing
#
#The pickle file contains a dictionary with structure:
#D=dict(system=dict(Porb=[],
#                   eccentricity=[],
#                   Wdisk=[],
#                   logQ=[],
#                   primary_mass=[],
#                   age=[],
#                   feh=[],
#                   Pspin=[])
#                   )

import numpy
import utils
import pickle

s=['1', '8', '12', '13', '17', '20', '25', '28', '32', '36', '39', '43', '44', '47', '48', '50', '54', '56', '57', '67', '70', '73', '76', '79', '80', '81', '83', '84', '85', '86', '88', '92', '93', '94', '95', '96', '106', '109', '120', '123', '126', '137']

s=['1']
params=dict()
parameters=['Porb','eccentricity','Wdisk','logQ','primary_mass','age','feh','Spin']


for system in s:
    D=dict()
    print(f'\nSystem={system}')
    delete_q=[]
    with open('deletion_checks/final_results.txt','r') as f:
        for lines in f:
            x=lines.split()
            s=x[0].split('_')[0]
            if s==system:
                if x[6]=='True':
                    delete_q.append(float(x[2]))
                else:delete_q.append(8.353237442958045)

    print(delete_q)

    logQ_values=utils._get_filled_chain(system)
    delete_idx=numpy.zeros(1)
    for q in delete_q:
        idx=numpy.where(logQ_values==q)
        delete_idx=numpy.concatenate((delete_idx,idx[0]),axis=None)
    delete_idx=delete_idx.astype(int)
    logQ=numpy.delete(logQ_values,delete_idx)
    for p in parameters:
        if p=='logQ':
            D[p]=logQ
        else:
            param_values=utils._get_filled_chain(system,p)
            param_values=numpy.delete(param_values,delete_idx)
            D[p]=param_values
    params[system]=D  
    
with open('complete_chains.pickle','wb') as f:
    pickle.dump(params,f)