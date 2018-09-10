class data_set :

    def __init__(self,**kwargs):

        for name,value in kwargs.items():
            setattr(self,name,value)



#rho_star = Structure(value=1.41,uncertainty=0.1)
data = data_set(x = data_set(a=1,b=2),
                y = data_set(a=3,b=4),
                z = data_set(a=5,b=6)
                )


print(data.z.a)