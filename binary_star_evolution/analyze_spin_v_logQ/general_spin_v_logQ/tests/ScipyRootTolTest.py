import scipy
from scipy import optimize

def fun(x):

    print('Trying x0=%s, x1=%s' %(repr(x[0]),repr(x[1])))

    x0 = x[0]  + 0.5 * (x[0] - x[1])**3 - 1.0
    x1 = 0.5 * (x[1] - x[0])**3 + x[1]

    print(x0,x1)

    return  x0,x1



sol = optimize.root(fun, [0, 0], method='hybr',options={'xtol':1e-3})
print(sol.x)
