import glob
import numpy

with open('period_dependence/new_stop_systems.txt','r') as f:
    converged_systems = f.read().split('\n')
converged_systems = [sk for sk in converged_systems if sk != '']

with open('aglib_error.txt','r') as f:
    error_systems = f.read().split('\n')
error_systems = [es.split('/')[0].split('_')[1] for es in error_systems]


all_systems = glob.glob('*.h5')
all_systems = [k.split('.')[0].split('_')[1] for k in all_systems]

remaining_error_systems = [k for k in error_systems if k not in converged_systems]
for i in range(0, len(remaining_error_systems), 8):
    print(' '.join(remaining_error_systems[i:i+8]))

print('\n')
remainder_good_systems = [g for g in all_systems if g not in converged_systems and g not in error_systems]
for i in range(0, len(remainder_good_systems), 8):
    print(' '.join(remainder_good_systems[i:i+8]))


