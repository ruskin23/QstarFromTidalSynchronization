import scipy as sc
import matplotlib.pyplot as plt

pi=sc.pi

class normal_distribution:

    #Creates a normal distrubution given a value for mu and sigma

    def __init__(self,mu,sigma):
        self.mu = mu
        self.sigma = sigma

    def distribution(self, value):

        A = 1/(self.sigma*sc.sqrt(2*pi))
        B = -0.5*((value - self.mu)/self.sigma)**2

        return A*sc.exp(B)



