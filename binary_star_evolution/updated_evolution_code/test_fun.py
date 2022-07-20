import numpy
from scipy import optimize

def test_fun(in_vals):

    print('testing values x = {} y = {}'.format(in_vals[0],in_vals[1]))
    
    a = numpy.sqrt(in_vals[0])
    b = numpy.sqrt(in_vals[1]) - 4
    
    # print('Calculated Values a = {} b = {}'.format(a,b))

    #return [numpy.sqrt(in_vals[0]), numpy.sqrt(in_vals[1]) - 4]
    print('Solution = {}'.format(numpy.sqrt(a**2 + b**2)))
    return  numpy.sqrt(a**2 + b**2)

initial_guess=[0.001,4]
# bounds=optimize.Bounds(numpy.array([0.0,1.0]),numpy.array([1.0,10.0]),keep_feasible=True)
sol = optimize.minimize(test_fun,
                    initial_guess,
                    # bounds=bounds,
                    method='Nelder-Mead'
                    )

print(sol.x)

