import numpy
from scipy import special
from scipy.stats import rv_continuous

def erf_fun(x):
    return 0.5*(1+special.erf(x/numpy.sqrt(2)))


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
        super().__init__(a=(self._samples.min()),
                        b=(self._samples.max()))
