import scipy as sc
import matplotlib.pyplot as plt

pi=sc.pi

class normal_distribution:

    #Creates a normal distribution given a value for mu and sigma

    def __init__(self,mu,sigma):
        self.mu = mu
        self.sigma = sigma

    def value_range(self):

        n=5
        min = self.mu - n*self.sigma
        max = self.mu + n*self.sigma
        v_array = sc.linspace(min, max, 1e+6)

        return v_array


    def distribution(self):


        A = 1/(self.sigma*sc.sqrt(2*pi))
        B = -0.5*(((self.value_range() - self.mu)/self.sigma)**2)

        return A*sc.exp(B)



def proposal_function(mu_value,sigma_prior,sigma_step):

    prior_function = normal_distribution(mu_value,sigma_prior).distribution()
    step_function = normal_distribution(mu_value,sigma_step).distribution()

    return prior_function*step_function

#picking a random value from proposal function

mu = 5668
sigma = 0.2

norm = normal_distribution(mu,sigma)

x = norm.value_range()
y = norm.distribution()


print (x)
print(y)

plt.plot(x,y)
plt.show()