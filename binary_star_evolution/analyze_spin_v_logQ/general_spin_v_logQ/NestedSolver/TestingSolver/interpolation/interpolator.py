import numpy
import scipy
import matplotlib.pyplot as plt
from scipy import interpolate


ei=[0.184,0.23800000000000002,0.24893609,0.265,0.29200000000000004]
pi=[16.625848635835336,16.940882141234365,17.015165351980404,17.129026858865828,17.408927806048062]

ef=[0.10031356866009071,0.1567115763857638,0.1672223320851414,0.18241389911336503,0.20292304359933191]

f=interpolate.interp2d(ei,pi,ef)

e_a=numpy.linspace(0.184,0.29200000000000004,10000)
p_a=numpy.linspace(16.625848635835336,17.408927806048062,10000)

ef_a=f(e_a,p_a)

print(ef_a)

ef_0=0.184

diff=ef_a[0,:]-ef_0

zero_crossing=numpy.where(numpy.diff(numpy.sign(diff)))[0]

print(zero_crossing)
print(e_a[zero_crossing])
print(p_a[zero_crossing])


"""
x1 = 0.184
x2 = 0.256
x3 = 0.328

f1 = 0.10031339436641426
f2 = 0.17382927564081121
f3 = 0.23368314775772145

x=numpy.linspace(x1,x3,10000)

q1 = (((x-x2)*(x-x3))/((x1-x2)*(x1-x3)))*f1
q2 = (((x-x1)*(x-x3))/((x2-x1)*(x2-x3)))*f2
q3 = (((x-x1)*(x-x2))/((x3-x1)*(x3-x2)))*f3

q=q1+q2+q3

q0=0.184

diff=q-q0

zero_crossing=numpy.where(numpy.diff(numpy.sign(diff)))[0]
print('Values:')
print(x[zero_crossing])
print(q[zero_crossing])
"""
