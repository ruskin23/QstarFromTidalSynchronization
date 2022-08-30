import numpy
import pickle
from scipy.stats import rv_continuous
import matplotlib.pyplot as plt
import corner

from utils import erf_fun
from utils import kernel_gauss

from pathlib import Path
from directories import directories

home_dir=str(Path.home())
path=directories(home_dir)


class kernel_gauss():

    def __init__(self,bandwidth):

        self.bandwidth=bandwidth


    def pdf(self,u):

        y = u/self.bandwidth

        K_h = ( 1
                /
                (numpy.sqrt(2*numpy.pi) * self.bandwidth) 
                ) * numpy.exp(-0.5*y*y)

        return K_h
    

    def cdf(self,u):

        K_c = erf_fun(u/self.bandwidth)

        return K_c 
        

class DiscreteSampling(rv_continuous):

    def _kernel_arg(self,x):
        rhs, lhs = numpy.meshgrid(self._samples, x)
        return numpy.squeeze(lhs - rhs)

    def _cdf(self,x):
        
        return numpy.average(self.kernel.cdf(self._kernel_arg(x)),
                            axis=-1,
                            weights=self.weights
                            )

    def _pdf(self,x):

        return numpy.average(self.kernel.pdf(self._kernel_arg(x)),
                            axis=-1,
                            weights=self.weight
                            )


    def get_norm(self,x):

        y = self._kernel_arg(x)/self._width

        Z = 1 / (self._width*numpy.sqrt(2*(numpy.pi)))

        return Z*numpy.exp(-0.5*y*y)


    def set_weights(self,weights):
        self.weights=weights

    def __init__(self,
                 samples,
                 width,
                 weights=None):

        self._samples=samples
        self._width=width
        self.weights=weights
        self.kernel=kernel_gauss(width)
        super().__init__(a=-numpy.inf,
                         b=numpy.inf)

def get_sample(_quantities,
                _samples,
                _bandwidths,
                unit_cube_iter
                ):

    weights=numpy.ones(len(_samples[_quantities[0]]))

    sampled_value=[]

    for Q in _quantities:

        DS=DiscreteSampling(_samples[Q],_bandwidths[Q])

        DS.set_weights(weights)

        sampled=DS.ppf(next(unit_cube_iter))

        weights*=DS.get_norm(sampled)

        sampled_value.append(sampled)
    
    return sampled_value

def get_system_params(systemKIC):

    system_chains=numpy.load(path.current_directory+'/catalog/samples/chains/'+systemKIC+'.npz')
    all_chains=numpy.transpose(system_chains['thinned_chain'])


    M_samples=all_chains[0]
    Q_samples=all_chains[1]
    Z_samples=all_chains[2]
    t_samples=all_chains[3]
    esinw_samples=all_chains[8]
    ecosw_samples=all_chains[9]
    
    e_samples=numpy.sqrt(ecosw_samples**2 + esinw_samples**2)

    _quantities=['Z','M','Q','t','e']

    _samples={'M':M_samples,
              'Q':Q_samples,
              'Z':Z_samples,
              't':t_samples,
              'e':e_samples}

    with open(path.current_directory+'/bandwidth.pickle','rb') as f:
        bw=pickle.load(f)
    bandwidths=bw[systemKIC]['ISJ']
    h_M=bandwidths[0]
    h_Q=bandwidths[1]
    h_Z=bandwidths[2]
    h_t=bandwidths[3]
    with open(path.current_directory+'/e_bandwidth.pickle','rb') as f:
        h_e=pickle.load(f)
    h_e=h_e[systemKIC]

    _bandwidths={'M':h_M,
                 'Q':h_Q,
                 'Z':h_Z,
                 't':h_t,
                 'e':h_e}

    return _quantities,_samples,_bandwidths

if __name__=='__main__':

    # systemKIC='9892471'
    
    with open('catalog/filtering/nominal_value_catalog_Iconv_cutoff.txt','r') as f:
        next(f)
        for lines in f:
            x=lines.split()
            systemKIC=x[1]
            print(systemKIC)
            
            _quantities,_samples,_bandwidths = get_system_params(systemKIC)

            sampled_values=[]

            for i in range(1):

                unit_cube=numpy.random.rand(5)
                unit_cube_iter=iter(unit_cube)

                sampled_values = numpy.append(sampled_values,get_sample(_quantities,_samples,_bandwidths,unit_cube_iter))
                
            print(sampled_values)
            sampled_values=numpy.reshape(sampled_values,(len(sampled_values)//5,5))

            figure=corner.corner(sampled_values)
            plt.savefig('compare_plots_bandwidths/system_{}/corner_ISJ_{}.png'.format(systemKIC,systemKIC))


