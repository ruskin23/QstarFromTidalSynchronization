import random
import scipy
from scipy.stats import norm



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


arr = []

for i in range(10):
    arr.append(i)

with open('test_3.txt','w') as file:
    for x in arr:
        file.write('%s ' %x)





print ("\n_______________________________________________________________________________________________________________")
print ("\n_______________________________________________________________________________________________________________")

if __name__ == '__main__':

    observation_data = dict(
                        age=dict(value=4.6, sigma=3.0),
                        teff=dict(value=5922.0, sigma=200.0),
                        feh=dict(value=-0.06, sigma=0.11),
                        semimajor=dict(value=0.05917, sigma=5.576e-09),
                        Porb=dict(value=5.2663825, sigma=3.7e-06),
                        Pdisk=dict(value=2*scipy.pi / 1.4, sigma=0.1)
                    )

    observed_Pspin = dict(
                        value=7.2,
                        sigma=0.1
                    )

    fixed_parameters = dict(
                        disk_dissipation_age=5e-3,
                        planet_formation_age=5e-3,
                        wind=True,
                        wind_saturation_frequency=2.54,
                        diff_rot_coupling_timescale=5e-3,
                        wind_strength=0.17,
                        inclination=scipy.pi/2

    )

    proposed_step = dict(
                        age_step=3.0,
                        teff_step=100.0,
                        feh_step=0.1,
                        semimajor_step=0.00001,
                        Porb_step=0.00001,
                        Pdisk_step=0.1,
                        logQ_step=0.5
                    )


    logQ = dict(
                min=4,
                max=6
            )

    prior = 1.0

    parameter_set = dict()

    for (name_obs, value_obs), (name_step, value_step) in zip(observation_data.items(),proposed_step.items()):
        print ("value =, step =  ", value_obs['value'], value_step)
        parameter_set[name_obs] = scipy.stats.norm.rvs(loc=value_obs['value'], scale=value_step)
        print ("parameter_values = ", parameter_set[name_obs])

    for (key_obs, value_obs), (key_parameter, value_parameter) in zip(observation_data.items(),parameter_set.items()):
        prior = scipy.stats.norm(value_obs['value'], value_obs['sigma']).pdf(value_parameter)
        print("obs_prior_name =, obs_prior_value = ",  key_obs, value_obs['value']   )
        print("par_prior_name =, par_prior_value = ",   key_parameter, value_parameter  )
        print ("prior = ", prior)



    propose_value = scipy.stats.norm.rvs(loc=0.05917, scale=1e-8)

    print (propose_value)

    check_prior = scipy.stats.norm(0.05917,5.576e-09).pdf(propose_value)
    print (check_prior)