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



x = test(observables,proposed_step)

x.ini()
x.disp()
