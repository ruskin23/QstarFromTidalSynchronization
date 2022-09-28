import numpy
from scipy import integrate
from scipy import special
import time

V=numpy.array([[ 9.80666324e+03, -4.88629693e+00, -8.44121472e-01],
               [ -4.88629693e+00,  7.51372889e-01, -4.23434050e-03],
               [ -8.44121472e-01, -4.23434050e-03, 6.23305007e-01]])


Vc=numpy.array([[-1.77480090e+01,  4.00299817e-01,  7.88742696e+00],
                [-1.01940598e+00, -1.76922521e-03,  1.80370153e-01],
                [-2.34457161e+00,  2.03731758e-01, -3.16692582e-01]])



def integration(theta_mean_vector):
    x0=theta_mean_vector[0]
    y0=theta_mean_vector[1]
    z0=theta_mean_vector[2]
    return (-0.029099*special.erf(6.69861 - 0.558217*z0)*special.erf(70.0808 - 70.0238*x0
                                                                  + 0.0348903*y0 -
                                                                  0.00602739*z0)*special.erf(0.23186
                                                                                             -
                                                                                             0.611938*y0
                                                                                             +
                                                                                             0.00311609*z0) +
         0.029099*special.erf(6.69861 - 0.558217 *z0)*special.erf(0.0569769 - 70.0238
                                                                  *x0 + 0.0348903 *y0 -
                                                                  0.00602739
                                                                  *z0)*special.erf(0.23186
                                                                                   -
                                                                                   0.611938
                                                                                   *y0 +
                                                                                   0.00311609
                                                                                   *z0) +
         0.029099*special.erf(2.79109 - 0.558217 *z0)*special.erf(70.0386 - 70.0238*x0 +
                                                                  0.0348903*y0 -
                                                                  0.00602739*z0)*special.erf(0.253672
                                                                                             -
                                                                                             0.611938*y0
                                                                                             +
                                                                                             0.00311609*z0) -
         0.029099*special.erf(2.79109 - 0.558217*z0)*special.erf(0.0147852 - 70.0238*x0
                                                                 + 0.0348903*y0 -
                                                                 0.00602739*z0)*special.erf(0.253672
                                                                                            -
                                                                                            0.611938*y0
                                                                                            +
                                                                                            0.00311609*z0) +
         0.029099*special.erf(6.69861 - 0.558217*z0)*special.erf(69.9398 - 70.0238*x0 +
                                                                 0.0348903*y0 -
                                                                 0.00602739*z0)*special.erf(2.70409
                                                                                            -
                                                                                            0.611938*y0
                                                                                            +
                                                                                            0.00311609*z0) -
         0.029099*special.erf(2.79109 - 0.558217*z0)*special.erf(69.8976 - 70.0238*x0 +
                                                                 0.0348903*y0 -
                                                                 0.00602739*z0)*special.erf(2.7259
                                                                                            -
                                                                                            0.611938*y0
                                                                                            +
                                                                                            0.00311609*z0) +
         0.029099*special.erf(6.69861 - 0.558217*z0)*special.erf(2.70409 - 0.611938*y0 +
                                                                 0.00311609*z0)*special.erf(0.0839797
                                                                                            +
                                                                                            70.0238
                                                                                            *x0
                                                                                            -
                                                                                            0.0348903*y0
                                                                                            +
                                                                                            0.00602739*z0) -
         0.029099*special.erf(2.79109 - 0.558217*z0)*special.erf(2.7259 - 0.611938*y0 +
                                                                 0.00311609*z0)*special.erf(0.126171
                                                                                            +
                                                                                            70.0238
                                                                                            *x0
                                                                                            -
                                                                                            0.0348903*y0
                                                                                            +
                                                                                            0.00602739*z0))





I=0
with open('../../../mcmc_mass_age/samples/updated_samples/MassAgeFehSamples_137.txt','r') as f:
    next(f)
    start_time=time.time()
    for lines in f:
        x=lines.split()
        phi_vector=numpy.array([float(x[k]) for k in range(3)]) - numpy.array([1.044975209491866,5.327761315229749,-0.10134764748374189])
        VV=numpy.matmul(numpy.linalg.inv(V).T,numpy.matmul(Vc,phi_vector))
        theta_mean_vector=numpy.array([0.011,4.1,8.442844284428443])-VV

        I=I+integration(theta_mean_vector)
        #func = lambda x1,x2,x3: numpy.exp(-0.5*numpy.matmul((numpy.array([x1,x2,x3])-theta_mean_vector).T,
        #                                                    numpy.matmul(V,
        #                                                                 (numpy.array([x1,x2,x3])-theta_mean_vector)
        #                                                                 )
        #                                                    )
        #                                  )
        #lim2=[0,1]
        #lim3=[2*numpy.pi/14,2*numpy.pi/1.4]
        #lim4=[5,12]
        #lim=[lim2,lim3,lim4]
        #I=I+integrate.nquad(func,lim)[0]
print(I)
print(time.time()-start_time)

"""
N=100

porb=numpy.linspace(7,9,N)
eccentricity=numpy.linspace(0,1,N)
Wdisk=numpy.linspace(2*numpy.pi/14,2*numpy.pi/1.4,N)
logQ=numpy.linspace(5,12,N)

mean_vector=numpy.array([8.20292763,0.00914629,4.03820967,8.39650154])


V=numpy.array([[ 1.22018405e+08,  4.07390464e+04,  6.02997425e+01,-8.41363730e+00],
               [ 4.07390464e+04,  9.80666324e+03, -4.88629693e+00, -8.44121472e-01],
               [ 6.02997425e+01, -4.88629693e+00,  7.51372889e-01, -4.23434050e-03],
               [-8.41363730e+00, -8.44121472e-01, -4.23434050e-03, 6.23305007e-01]])

I=[]


#for p in porb:

for e in eccentricity:
    for w in Wdisk:
        for q in logQ:
            theta_vector=numpy.array([e,w,q])-mean_vector
            I=numpy.append(I,numpy.exp(-0.5*numpy.matmul(theta_vector.T,numpy.matmul(V,theta_vector))))

print(max(I))
print(min(I))
"""
