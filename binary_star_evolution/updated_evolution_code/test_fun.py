import numpy
from scipy import optimize

def test_fun(in_vals):

    print('testing values x = {} y = {}'.format(in_vals[0],in_vals[1]))
    
    a = numpy.sqrt(in_vals[0])
    b = numpy.sqrt(in_vals[1]) - 4
    
    print('Calculated Values a = {} b = {}'.format(a,b))

    if numpy.isnan(a) or numpy.isnan(b):
        return numpy.array([numpy.inf,numpy.inf])

    return [numpy.sqrt(in_vals[0]), numpy.sqrt(in_vals[1]) - 4]

initial_guess=[-0.1,2.3]
sol = optimize.least_squares(test_fun,
                    initial_guess
                    )

print(sol.x)

