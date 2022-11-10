import glob
import json
import numpy

converged_list_kalo_plots = ['11147276', '10091257', '5022440', '5802470', '4815612', '3241344', '10031409', '7838639', '8957954', '4285087', '12004679', '9971475', '10385682', '7362852', '8543278', '6927629', '8364119', '6949550', '3348093', '4839180', '5652260', '11232745', '8984706', '5393558', '9353182', '10936427', '4678171', '7987749', '10330495', '10711913', '10215422', '4773155']

all_systems = glob.glob('*.h5')
all_systems = [k.split('.')[0].split('_')[1] for k in all_systems]


def kalo_converged_list():

    not_converged = []
    converged = []
    for kic in all_systems:
        if kic not in converged_list_kalo_plots:
            not_converged.append(kic)
        else: converged.append(kic)

    print('\nNOT CONVERGED:')
    for i in range(0, len(not_converged), 8):
        print(' '.join(not_converged[i:i+8]))
    
    print('\nCONVERGED:')
    for i in range(0, len(converged), 8):
        print(' '.join(converged[i:i+8]))



def local_converged_list():

    not_converged = []
    converged = []
    covnergence_dict = dict()
    for kic in all_systems:
        burnins = numpy.ones((50,4), dtype=int)
        with open(f'period_dependence/{kic}.txt', 'r') as f:
            next(f)
            for i, lines in enumerate(f):
                x = lines.split()
                burnins[i] = numpy.array([int(x[k].split('/')[0]) for k in [3, 6, 9, 12]])
                max_step = int(x[3].split('/')[1])
                if 'nan' in x and kic not in not_converged: not_converged.append(kic)
        if kic not in not_converged:
            if max_step - numpy.max(burnins)<100: not_converged.append(kic)
            else: converged.append(kic)

        covnergence_dict[kic] = dict()
        covnergence_dict[kic]['max_burn_in'] = repr(numpy.max(burnins))
        covnergence_dict[kic]['max_step'] = repr(max_step)
        covnergence_dict[kic]['step_diff'] = repr(max_step - numpy.max(burnins))
        if kic in converged: covnergence_dict[kic]['converged'] = 'True'
        else: covnergence_dict[kic]['converged'] = 'False'


    with open("convergence.json", "w") as outfile:
        json.dump(covnergence_dict, outfile, indent=4)

    # for i in range(0, len(not_converged), 8):
    #     print(' '.join(not_converged[i:i+8]))

    for i in range(0, len(not_converged), 8):
        print(' '.join(not_converged[i:i+8]))


if __name__ == '__main__':

    kalo_converged_list()

# 5022440 4678171 8957954 7838639 4815612 7987749 7362852 5393558
# 10936427 4839180