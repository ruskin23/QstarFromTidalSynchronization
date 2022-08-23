import errno
import numpy
from scipy import optimize

# def test_fun(in_vals):

#     print('testing values x = {} y = {}'.format(in_vals[0],in_vals[1]))
    
#     a = numpy.sqrt(in_vals[0])
#     b = numpy.sqrt(in_vals[1]) - 4
    
#     # print('Calculated Values a = {} b = {}'.format(a,b))

#     #return [numpy.sqrt(in_vals[0]), numpy.sqrt(in_vals[1]) - 4]
#     print('Solution = {}'.format(numpy.sqrt(a**2 + b**2)))
#     return  numpy.sqrt(a**2 + b**2)


# bounds=optimize.Bounds(numpy.array([0.0,1.0]),numpy.array([1.0,10.0]))
# initial_simplex=numpy.array([[1,2],[4,2],[4,7]])
# initial_guess = (numpy.array([1,2])+numpy.array([4,2])+numpy.array([4,7]))/3
# sol = optimize.minimize(test_fun,
#                     initial_guess,
#                     method='Nelder-Mead',
#                     bounds=bounds,
#                     options={'initial_simplex' : initial_simplex}
#                     )

# print(sol.x)



# class Num:
#     def __init__(self,num,n2=4):
#         self.n1 = num
#         self.n2=n2
#     def square(self):
#         print(self.n2**2)

# class Num2(Num):
#     def show(self):
#         print(self.n1)

# N=Num2(2,n2=9)
# N.square()


b=34.3
atol=1e-5
rtol=1e-6

err=atol + (rtol * abs(b))



cache=dict()


cache[(1,2)]={'a':1.2,'b':4.2}
cache[(3,4)]={'a':1.2,'b':4.2}


print(list(cache.keys())[-1])

name='a'
if (1,2) in cache:
    print(cache[(1,2)]['a'])

