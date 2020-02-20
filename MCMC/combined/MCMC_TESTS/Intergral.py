import numpy
from scipy import special
from scipy import integrate
import matplotlib.pyplot as plt
import itertools


logQ_mean=7.125744153362705
t_mean=4.049835253026934


sigma_logQ=0.2
sigma_t=0.8


logQ_min=5.0
logQ_max=12.0


def intergand(t,logQ):

    I = numpy.exp( - ((logQ-logQ_mean)/(sigma_logQ))**2 - (((t*logQ)-(t_mean*logQ_mean))/(sigma_t*sigma_logQ))**2 )

    return I

def EvaluateIntegral(t):
    A = (1/sigma_logQ**2)*(1 + (t**2)/(sigma_t**2))
    B = 2*(1 + (t*t_mean)/(sigma_t**2))*((logQ_mean)/(sigma_logQ**2))
    C = (1 + ((t_mean**2)/((sigma_t**2))))*((logQ_mean**2)/(sigma_logQ**2))

    alpha = (C/A) - (B**2/A**2)

    Omega = numpy.exp(A*alpha)*numpy.sqrt(numpy.pi/A)/2

    x1 = numpy.log(10)*logQ_min
    x2 = numpy.log(10)*logQ_max

    y1 = x1 - (B/(2*A))
    y2 = x2 - (B/(2*A))

    z1 = numpy.sqrt(A)*y1
    z2 = numpy.sqrt(A)*y2


    I = (numpy.sqrt(numpy.pi)/2)*Omega*(special.erf(z2)-special.erf(z1))
    """

    B = t/(2*sigma_t)
    C = ((t-t_mean)/sigma_t)*(logQ_mean/sigma_logQ)

    alpha = -B**2 + C

    Omega = (numpy.sqrt(numpy.pi)/2)*numpy.exp(-alpha)*sigma_logQ

    x1 = numpy.log(10)*logQ_min
    x2 = numpy.log(10)*logQ_max

    y1 = (x1 - logQ_mean)/sigma_logQ
    y2 = (x2 - logQ_mean)/sigma_logQ

    z1 = y1 + B
    z2 = y2 + B

    I = Omega*(special.erf(z2)-special.erf(z1))
    """


    print('age = ',t)
    print('A = ',A)
    print('B = ',B)
    print('C = ',C)
    print('alpha = ',alpha)
    print('Omega = ',Omega)
    print('z1 = ',z1)
    print('z2 = ',z2)
    print('erf1 = ',special.erf(z1))
    print('erf2 = ',special.erf(z2))
    print('\n')

    return I


def EvaluateProbability(t,m):

    I = lambda logQ: intergand(t,logQ)
    print('for t = {} I = {}'.format(t,integrate.quad(I,logQ_min,logQ_max)[0]))
    P =    m*numpy.exp(-((t-t_mean)**2)/(sigma_t**2))*integrate.quad(I,logQ_min,logQ_max)[0]

    return P


def cummulative_distribution(x):

    a=[]
    v=0
    value=[]

    for tuples in x:
        a=numpy.append(a,tuples[0])
        v=v+tuples[1]
        value=numpy.append(value,v)

    return list(zip(a,value/max(value)))



def MarginalizeLogQ(age,multiplicity):

    N=0
    modified_multiplicity=[]
    for a,m in zip(age,multiplicity):
        P=EvaluateProbability(a,m)
        modified_multiplicity=numpy.append(modified_multiplicity,P)
        N=N+P

    modified_multiplicity=modified_multiplicity/N

    age_tuple=list(zip(age,modified_multiplicity))
    age_tuple=sorted(age_tuple, key=lambda tup: tup[0])

    age_cummulative=cummulative_distribution(age_tuple)

    return age_cummulative

def MarginalizeAge(age,multiplicity,logQ):

    s=0

    for t,m in zip(age,multiplicity):

        m_t = multiplicity
        arg = ((t-t_mean)/sigma_t)**2 + ((t*logQ - t_mean*logQ_mean)/(sigma_t*sigma_logQ))**2

        s = s + m*numpy.exp(-arg)


    return s

def CompareLogQ(age,multiplicity):

    logQ_uniform=numpy.linspace(logQ_min,logQ_max,1000)
    logQ_dist=numpy.exp(-((logQ_uniform-logQ_mean)/sigma_logQ)**2)*MarginalizeAge(age,multiplicity,logQ_uniform)

    logQ_tuple=list(zip(logQ_uniform,logQ_dist))

    logQ_cummulative=cummulative_distribution(logQ_tuple)


    plt.scatter(*zip(*logQ_cummulative),color='r')



    mcmc_file='/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/combined/MCMC_TESTS/TEST5/AccetedParameters.txt'

    logQ_mcmc=[]
    with open(mcmc_file,'r') as f:
        next(f)
        for lines in f:
            x=lines.split()
            logQ_mcmc=numpy.append(logQ_mcmc,float(x[4]))

    logQ_mcmc_tuple=[(x, len(list(y))) for x, y in itertools.groupby(logQ_mcmc)]
    logQ_mcmc_tuple=sorted(logQ_mcmc_tuple, key=lambda tup: tup[0])

    logQ_mcmc_cummulative=cummulative_distribution(logQ_mcmc_tuple)

    plt.plot(*zip(*logQ_mcmc_cummulative))
    plt.show()


if __name__=='__main__':

    age=[]
    multiplicity=[]
    N=0

    filename='/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/mcmc_mass_age/samples/updated_samples/MassAgeFehSamples_54.txt'

    with open(filename,'r') as f:
        next(f)
        for lines in f:
            x=lines.split()
            age=numpy.append(age,float(x[1]))
            multiplicity=numpy.append(multiplicity,float(x[3]))


    age_cummulative=MarginalizeLogQ(age,multiplicity)

    age_cummulative=cummulative_distribution(age_cummulative)
    plt.scatter(*zip(*age_cummulative))


    mcmc_file='/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/combined/MCMC_TESTS/TEST5/AccetedParameters.txt'

    age_mcmc=[]
    with open(mcmc_file,'r') as f:
        next(f)
        for lines in f:
            x=lines.split()
            age_mcmc=numpy.append(age_mcmc,float(x[6]))

    age_mcmc_tuple=sorted([(x, len(list(y))) for x, y in itertools.groupby(age_mcmc)], key=lambda tup: tup[0])

    age_mcmc_cummulative=cummulative_distribution(age_mcmc_tuple)

    plt.plot(*zip(*age_mcmc_cummulative))

    plt.show()


    plt.show()

    """
    t_value=numpy.linspace(3.0,5.0,10)
    colors=['r','g','b','c','m','y','k','lime','gray','gold']
    for t,c in zip(t_value,colors):

        logQ_uniform=numpy.linspace(logQ_min,logQ_max,1000)
        #integrand = numpy.exp( - ((logQ_uniform-logQ_mean)/(sigma_logQ))**2 - (((t*logQ_uniform)-(t_mean*logQ_mean))/(sigma_t*sigma_logQ))**2 )

        integrand=reduced_integrand(t)
        plt.semilogy(logQ_uniform,integrand,color=c,label=t)

        plt.plot()
    plt.ylim(10e-10,1)
    plt.legend()
    plt.show()
    """
