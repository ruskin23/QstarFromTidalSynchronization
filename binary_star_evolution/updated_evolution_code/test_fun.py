import numpy
from scipy import optimize

def test_fun(in_vals):

    print('testing values x = {} y = {}'.format(in_vals[0],in_vals[1]))


    return numpy.sqrt(in_vals[0]), numpy.sqrt(in_vals[1]) - 4

initial_guess=[-0.1,2.3]
sol = optimize.root(test_fun,
                    initial_guess,
                    method='lm',
                    tol=1e-6,
                    options={'maxiter':20}
                    )

print(sol.x)

