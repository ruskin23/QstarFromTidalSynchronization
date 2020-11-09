import matplotlib.pyplot as plt
import numpy

emin=numpy.linspace(0,0.42,100)

N=[16,20,25,30,45]


for n in N:

    delta_e=(0.43-emin)/n

    plt.plot(emin,delta_e,label=str(n))
plt.legend()
plt.show()
