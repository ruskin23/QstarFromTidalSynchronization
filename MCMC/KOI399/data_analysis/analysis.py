import pickle
import matplorlib.pyplot as plt
import numpy

def plot_output(fname,parameters,*argv):

    with open(fname,'r') as f:
        reader =csv.reader(f, dialect='excel-tab')
        for line in reader:
            param['iterations'].append(int(line[0]))
            for index,key in enumerate(parameters):
                param[key].append(float(line[index+1]))


    plt.scatter(param[argv],logQ,marker = '.')
    plt.show()




if __name__ == '__main__':

    parameters = [ 'age', 'teff', 'feh', 'Wdisk', 'logQ']

    #combine all formatted accepted parameter files
    in_fnames = ['formatted_accepted_parameters_'+repr(k)+'.txt' for k in range(2,8)]
    out_fname ='combined_accepted_parameters.txt'
    with open(out_fname,'w+') as fout:
        for name in in_fnames:
            with open(name,'r') as fin:
                for line in fin:
                    fout.write(line)


    #plot_parameters
    plot_output(out_fname,parameters,'age')









