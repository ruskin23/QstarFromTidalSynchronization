import matplotlib.pyplot as plt
import numpy
from scipy import stats
from utils import cummulative_distribution

clusters=['ganymede','stampede']

def _get_filename(system,cluster,instance):
#returns the filename with directory of accepted parameter file

    return '../'+cluster+'/MCMC_'+system+'/accepted_parameters_'+str(instance)+'.txt'

def _write_on_file(line,filename,option):
#write line on given filename
#option:'w' to write new file
#       'a' to append

    with open(filename,option) as f:
        f.write(line)

def _fill_parameters(chain_filename):
#creates a file filled_parameters.txt which repeats the missing parameters
#in accepted parameter file

    filled_filename='filled_parameters.txt'
    with open(chain_filename,'r') as f:
        next(f)
        counter=1
        for i,params in enumerate(f):
            x=params.split()
            if i==0:
                if x[0]!='1':raise ValueError
                else:
                    _write_on_file(params,filled_filename,'w')
                    saved_state=x[1:-1]
                    counter=counter+1
                    continue
            else:
                iteration=int(x[0])
                if iteration==counter:
                    _write_on_file(params,filled_filename,'a')
                    saved_state=x[1:-1]
                    counter=counter+1
                    continue
                else:
                    difference=iteration-counter
                    for _ in range(difference):
                        line='\t'.join([str(counter)]+saved_state+['\n'])
                        _write_on_file(line,filled_filename,'a')
                        counter=counter+1
                    _write_on_file(params,filled_filename,'a')
                    saved_state=x[1:-1]
                    counter=iteration+1


    return filled_filename

def get_cdf(parameter,filled_filename):
#returns cdf for the specified paremter in file of accepted parameters
#the filed_filename is the file where repeated steps is accounted for

    samples=[]
    with open(filled_filename,'r') as f:
        for lines in f:
            x=lines.split()
            samples=numpy.append(samples,float(x[4]))

    C=cummulative_distribution(samples)
    return C()

def get_chain(chain_filename):

    CHAIN=[]
    with open(chain_filename,'r') as f:
        next(f)
        for lines in f:
            x=lines.split()
            CHAIN=numpy.append(CHAIN,float(x[4]))
    return CHAIN

def get_stats(chain1,chain2):

    return stats.ks_2samp(chain1,chain2)

def break_chain(system,percent_break):

    complete_chain=[]
    for c in clusters:
        for i in range(5):
            _fill_parameters(_get_filename(system,c,i+1))
            with open('filled_parameters.txt','r') as f:
                for lines in f:
                    x=lines.split()
                    complete_chain=numpy.append(complete_chain,float(x[4]))

    split_chain=numpy.split(complete_chain,[int(len(complete_chain)*(percent_break/100))])

    print(get_stats(split_chain[0],split_chain[1]))


def all_chains(system):

    chain_filenames=[]
    for c in clusters:
        for i in range(5):
            name=_get_filename(system,c,i+1)
            chain_filenames.append(name)

    #Totally not optimal. definitely need more editing. But lengths are not that big so this is running fast for now
    for i in range(10):
        if i<5:i_filename='ganymede chain '+str(i+1)
        else:i_filename='stampede chain '+str(i-4)
        for j in range(i+1,10):
            if j<5:j_filename='ganymede chain '+str(j+1)
            else:j_filename='stampede chain '+str(j-4)
            comparison_between='For '+ i_filename + ' and ' + j_filename
            chain1=get_chain(chain_filenames[i])
            chain2=get_chain(chain_filenames[j])
            print(comparison_between,get_stats(chain1,chain2))

def plot_chains(parameter,chain_filenames):
#plot cdf of two different chains


    for i in range(2):
        filled_filename=_fill_parameters(chain_filenames[i])
        cdf=get_cdf(parameter,filled_filename)
        plt.plot(*zip(*cdf))
    plt.show()

if __name__ == '__main__':

    system='84'
    parameters=['Porb','eccentricity','Wdisk','logQ','primary_mass','age','feh','Pspin']

    #all_chains(system)
    #break_chain(system,70)
    chain1=_get_filename(system,'ganymede',5)
    chain2=_get_filename(system,'stampede',3)
    plot_chains('logQ',[chain1,chain2])
