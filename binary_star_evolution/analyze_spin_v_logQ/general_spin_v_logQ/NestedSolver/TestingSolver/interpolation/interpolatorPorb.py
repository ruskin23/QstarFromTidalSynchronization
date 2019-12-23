import numpy



x1=15.927
x2=16.927

y1=15.175577498495041
y2=16.248028035513023

m=(y2-y1)/(x2-x1)
c=(x2*y1-x1*y2)/(x2-x1)


x=numpy.linspace(x1,x2,1000000)
y=m*x+c

y0=15.927

diff=y-y0

zero_crossing=numpy.where(numpy.diff(numpy.sign(diff)))[0]

print(x[zero_crossing])
print(y[zero_crossing])

