import numpy


class get_arguments(object):

    def __init__(self,parameters):
        for key,item in parameters.items():
            setattr(self,key,get_arguments(item) if isinstance(item,dict) else item)

#parameters = {'primary_mass' : {'min' : 0.5, 'max' : 1.4}}


parameters = dict(
        primary_mass=dict(min=0.4,max=1.4),
        feh=dict(value=0.1,sigma=0.1),
        age=dict(min=1,max=13)
    )


x=get_arguments(parameters)

x.primary_mass.min = 4
