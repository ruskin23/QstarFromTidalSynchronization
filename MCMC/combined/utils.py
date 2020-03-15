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


