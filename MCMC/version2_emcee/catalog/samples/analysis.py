import numpy
import corner
import matplotlib.pyplot as plt
import sys


def kernel_gauss(x,x_i,h):

    y=(x-x_i)/h
    return (1/(numpy.sqrt(2*numpy.pi)))*numpy.exp(-(y**2)/2) 

def kde(x,x_array,h=None):

    n=len(x_array)
    if h is None:h=0.001#3.5*numpy.power(n,-1.0/3)

    f=0
    for i in range(n):
        f=f+kernel_gauss(x,x_array[i],h)
    return f/(n*h)


def get_eccentricity_col(all_chains):

    all_chains=numpy.transpose(all_chains)

    esinw_samples=all_chains[8]
    ecosw_samples=all_chains[9]

    e_samples=[]
    for i in range(len(esinw_samples)):
        f=esinw_samples[i]**2 + ecosw_samples[i]**2
        f=numpy.sqrt(f)
        e_samples.append(f)
    e_samples=numpy.array(e_samples)

    return e_samples


def eccentricity_plots(all_chains):

    sliced_chain=all_chains[:,8:10]
    figure=corner.corner(sliced_chain)
    plt.show()

def stellar_paramters_plots(all_chains):
    
    
    sliced_chain=all_chains[:,0:4]

    sliced_chain[:,[2,0]]=sliced_chain[:,[0,2]]
    sliced_chain[:,[2,1]]=sliced_chain[:,[1,2]]

    eccentricity_col=get_eccentricity_col(all_chains)
    eccentricity_col=numpy.reshape(eccentricity_col,[len(eccentricity_col),1])
    print(eccentricity_col.shape)
    sliced_chain=numpy.hstack((sliced_chain,eccentricity_col))

    print(numpy.shape(sliced_chain))

    figure=corner.corner(sliced_chain)
    # plt.show()
    plt.savefig(f'corner_{system}.png')



if __name__=='__main__':

    system = sys.argv[1]

    # A=numpy.load('chains/10198109.npz')
    A=numpy.load('chains/'+system+'.npz')

    all_chains=A['thinned_chain']

    stellar_paramters_plots(all_chains)

    params=['Msum','Q','z','age','d','E','Porb','tPE','esinw','ecosw','b','q11','q12','q21','q22','lnsigmaLC','lnsigmaSED','lnsigmaE','lnsigmad']

    # eccentricity_plots(all_chains)

    # all_chains=numpy.transpose(all_chains)

    # parameters=dict()
    # for i,param in enumerate(params):
    #     parameters[param]=all_chains[i]

    # #sum of masses
    # x=numpy.linspace(0.001,0.06,1000)
    # pdf_M=kde(x,parameters['z'])
    # plt.plot(x,pdf_M)
    # plt.show()

#make sure not breaking detailed balance for the step function: will not work
#integrate over parameters and convolve kernel to make continuous distribution: make sure not breaking aasumption of emcee like 1. 

