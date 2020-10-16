import random
import scipy
from scipy.stats import norm
import numpy
import matplotlib.pyplot as plt

class struct :

    def __init__(self,**kwargs):

        for name,value in kwargs.items():
            setattr(self,name,value)



class test:

    def disp(self):

        print (self.set)

    def ini(self):

        for i in range(10):

            self.set.append(i)


    def __init__(self,observables,proposed_step,set = None):


        self.observables = observables
        self.proposed_step = proposed_step
        if set is None:

            self.set = []




    observation_data = struct(
                        Teff = struct( value = 1, sigma = 0.1 ),
                        feh = struct(value=1, sigma=0.1),
                        rvk = struct(value=1, sigma=0.1),
                        inclination = struct(value=1, sigma=0.1)
    )

    fixed_parameters = struct(
                                disk_dissipatoin_age = 5e-3
    )

    proposed_step = struct(
                        Teff_proposed_step = 0.2,
                        feh_proposed_step = 0.2,
                        rvk_proposed_step = 0.3,
                        inclination_step = 0.2
    )





class test2:

    def __init__(self,p):

        for name,value in p.items():
            setattr(self,name,value)

    def printing(self):


        print (self.b['s'])



d = dict(a = dict(v = 1, s = 0.1),
         b = dict(v = 2, s = 0.2),
         c = dict(v = 3, s = 0.3))



p = dict(psa = 0.2,
         psb = 0.3,
         psc = 0.4
         )

instance = test2(d)
#instance.printing()
n = dict()

check = p

observation_data = dict(
    Teff=dict(value=1, sigma=0.1),
    feh=dict(value=6, sigma=0.2),
    rvk=dict(value=4, sigma=0.3),
    inclination=dict(value=23, sigma=0.1)
)

#for key,value in observation_data.items():

    #print (key,value['sigma'])

#for (kd,vd),(kp,vp) in zip(d.items(),p.items()):
    #print (kd)
    #print (vp)
#    n[ kd] = vp*4

#for x in check:
#        print (x)#print (d[x]['v'])

width = scipy.linspace(0,1,100)

dist = norm.pdf(width,loc=0.5,scale=0.1)

class test2:

    def __init__(self,x,y):

        self.a = x
        self.b = y

    def __call__(self):

        z = self.a*self.b

        return z



mean = 0
sigma = 1

#s = numpy.random.normal(loc = mean, scale = sigma, size = 1000)
s = scipy.stats.norm.rvs(loc = mean, scale = sigma, size = 100000)
plt.hist(s, bins = 30, normed=True)
plt.show()

#for i in range(1000):
#    s[i] = numpy.random.normal(loc = mean, scale = sigma)
#    print(s[i])



