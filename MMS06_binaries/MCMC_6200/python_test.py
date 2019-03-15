from scipy import optimize

def fun(x):

    print('in fun')
    print(x[0])
    print(x[1])
    c=0.5
    return [x[0]  + c * (x[0] - x[1])**3 - 1.0,
            c * (x[1] - x[0])**3 + x[1]]


sol = optimize.root(fun, [0, 0], method='hybr')
print(sol.x)

