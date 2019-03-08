import scipy
import numpy
import matplotlib.pyplot as plt




def fun(m):
    return (m**3)/((0.7+m)**2)


m = numpy.arange(0.01,0.7,0.001)

ms = fun(m)
c = [0.00032653311042704614]*len(ms)
print(ms - c)
plt.plot(m,ms)
plt.plot(m,c)
plt.show()
