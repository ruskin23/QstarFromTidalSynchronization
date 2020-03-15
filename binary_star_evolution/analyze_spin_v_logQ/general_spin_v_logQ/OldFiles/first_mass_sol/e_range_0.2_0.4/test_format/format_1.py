import os

directory='/home/ruskin/projects/QstarFromTidalSynchronization/binary_star_evolution/analyze_spin_v_logQ/general_spin_v_logQ/first_mass_sol/e_range_0.2_0.4/test_format'

for filename in os.listdir(directory):
    if filename.endswith('.txt'):

        a =  dict()
        with open(filename,'r') as f:
            for i,line in enumerate(f):
                a[i] = line

        with open(filename,'w') as f:
            f.write(a[0] + a[1] + a[2] + a[7] + a[8] + a[3] + a[4] + a[5] + a[6])

