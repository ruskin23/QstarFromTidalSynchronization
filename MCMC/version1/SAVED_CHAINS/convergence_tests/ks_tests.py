import matplotlib.pyplot as plt
import numpy
from scipy import stats
from utils import _cummulative_distribution
from utils import _get_filename
from utils import _write_on_file
from utils import _fill_parameters

clusters=['ganymede','stampede']


def get_cdf(filled_filename):
#returns cdf for the specified paremter in file of accepted parameters
#the filed_filename is the file where repeated steps is accounted for

    samples=[]
    with open(filled_filename,'r') as f:
        for lines in f:
            x=lines.split()
            samples=numpy.append(samples,float(x[4]))

    return _cummulative_distribution(samples)


def get_chain(chain_filename):
    CHAIN=[]
    with open(chain_filename,'r') as f:
        next(f)
        for lines in f:
            x=lines.split()
            CHAIN=numpy.append(CHAIN,float(x[4]))
    return CHAIN


def get_stats(chain1,chain2):
#returns ks 2 sample statistics using scipy ks_2samp function
#accepts 2 numpy.ndarray of chains

    return stats.ks_2samp(chain1,chain2)


def break_chain(system,percent_break):
    
    chain1=[]
    chain2=[]
    for c in clusters:
        for i in range(5):
            _fill_parameters(_get_filename(system,c,i+1))
            complete_chain=get_chain('filled_parameters.txt')
            for k,v in enumerate(complete_chain):
                if k<int(len(complete_chain)*(percent_break/100)): chain1=numpy.append(chain1,v)
                else:chain2=numpy.append(chain2,v)
            #split_chain=numpy.split(complete_chain,[int(len(complete_chain)*(percent_break/100))])
    print(get_stats(chain1,chain2))
    return [chain1,chain2]


def all_chains(system):
#Useless for now

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


def plot_chain(chain,label=None):
#plot cdf of two different chains

    cdf=_cummulative_distribution(chain)
    plt.plot(*zip(*cdf),label=label)


if __name__ == '__main__':

    system='43'

    #all_chains(system)
    
    chains=break_chain(system,50)
    for i,c in enumerate(chains):
        plot_chain(c,label='part'+str(i+1))
    complete_chain=numpy.concatenate([chains[0],chains[1]])
    plot_chain(complete_chain,label='complete')
    plt.legend()
    plt.show()    