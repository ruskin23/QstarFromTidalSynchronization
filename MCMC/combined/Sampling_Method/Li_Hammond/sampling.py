from utils import cummulative_distribution
import numpy
from scipy import special
import matplotlib.pyplot as plt
from scipy import integrate


class sampling_algo:


    def __init__(self,
                 parameters,
                 sample_filename,
                 rho_v):

        self.parameters=parameters
        self.sample_filename=sample_filename
        self.rho_v=rho_v


    def get_CDF(self):

        mass_samples=[]
        feh_samples=[]
        multiplicity=[]
        with open(sample_filename,'r') as f:
            next(f)
            for lines in f:
                x=lines.split()
                mass_samples=numpy.append(mass_samples,float(x[0]))
                feh_samples=numpy.append(feh_samples,float(x[2]))
                multiplicity=numpy.append(multiplicity,float(x[3]))

        self.CDF_mass=cummulative_distribution(mass_samples,multiplicity)()
        self.CDF_feh=cummulative_distribution(feh_samples,multiplicity)()



    def zi(self,
           v,
           key,
           error):

        F=(v-self.parameters[key]['value'])/self.parameters[key]['sigma']
        erf_value=0.5*(1+special.erf(F/numpy.sqrt(2)))
        if key=='mass':CDF=self.CDF_mass
        if key=='feh':CDF=self.CDF_feh
        CDF=numpy.array(CDF)
        C=CDF[:,1]-erf_value
        zero_crossings = numpy.where(numpy.diff(numpy.sign(C)))[0]

        print('for {}, v = {}, ZERO = {}'.format(key,v,zero_crossings))
        if not zero_crossings:
            if v<1:
                print('Setting zero_crossings for {} = 0'.format(key))
                zero_crossings=0
                print('for {}, v = {}, ZERO = {}'.format(key,v,zero_crossings))
            else: print('No SOlutions Found for key = {} because z={}'.format(key,zero_crossings))

        return CDF[zero_crossings,0]

    def integrand(self,
                  v1,
                  v2):

        print('\nv1 = {}, v2 = {}'.format(v1,v2))
        z1=self.zi(v1,'mass',1e-5)
        z2=self.zi(v2,'feh',1e-5)
        print('z1 = ',z1)
        print('z2 = ',z2)

        f_v=numpy.exp(-(v1**2 - 2*v1*v2 + v2**2)/(2*(1-self.rho_v**2)))
        f_v=(1.0/(2*numpy.pi*numpy.sqrt(1-self.rho_v**2)))*f_v
        print('f_v = ',f_v)
        print('I =', z1*z2*f_v)

        return z1*z2*f_v

    def integral(self):

        self.get_CDF()
        f = lambda v1,v2 :self.integrand(v1,v2)
        return integrate.dblquad(f,-1,0.5,lambda x:0.4, lambda x:1.2)


if __name__=='__main__':


    sample_filename='../../mcmc_mass_age/samples/updated_samples/MassAgeFehSamples_10.txt'
    rho_v=0.2

    parameters=dict(mass=dict(value=1.0,
                              sigma=0.5),
                    feh=dict(value=0.0,
                             sigma=0.3)
                    )

    I=sampling_algo(parameters,sample_filename,rho_v)
    print(I.integral())
    """
    mass=[]
    multiplicity=[]
    with open(sample_filename,'r') as f:
        next(f)
        for lines in f:
            x=lines.split()
            mass=numpy.append(mass,float(x[0]))
            multiplicity=numpy.append(multiplicity,float(x[3]))



    CDF=cummulative_distribution(mass,multiplicity)()
    v=numpy.linspace(-10,10,len(CDF))
    erf_value=special.erf(v)

    plt.plot(*zip(*CDF))
    plt.plot(v,erf_value)
    plt.show()
    """
