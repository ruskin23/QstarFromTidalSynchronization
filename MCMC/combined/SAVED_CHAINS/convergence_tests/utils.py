import numpy
import itertools
import matplotlib.pyplot as plt

class cummulative_distribution:

    def __init__(self,
                 samples,
                 mulitplicity=None,
                 plot=None):

        self.samples=samples
        self.mulitplicity=mulitplicity

    def sort_tuples(self):

        if self.mulitplicity is None: return sorted([(x, len(list(y))) for x, y in itertools.groupby(self.samples)], key=lambda tup: tup[0])
        else: return sorted(zip(self.samples,self.mulitplicity), key=lambda tup: tup[0])


    def __call__(self):


        sorted_samples=self.sort_tuples()

        sample_values=[]
        probability=[]
        p_sum=0.0

        for s in sorted_samples:
            sample_values=numpy.append(sample_values,s[0])
            p_sum=p_sum+s[1]
            probability=numpy.append(probability,p_sum)

        return list(zip(sample_values,probability/max(probability)))



class multivariate_gaussian:

    def multi(self,
              x_vector,
              x_mean,
              y_vector,
              y_mean,
              sigma_xy):

        arg=numpy.matmul(numpy.transpose(x_vector-x_mean),numpy.matmul(sigma_xy,(y_vector-y_mean)))
        return numpy.exp(-0.5*arg)


class  Norm:

    def N(self,
          value,
          loc=0.0,
          sigma=1.0):

        arg=(value-loc)/sigma
        return numpy.exp(-(arg**2)/2)

