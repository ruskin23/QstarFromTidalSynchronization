import numpy as np
import matplotlib.pyplot as plt
import corner


with open('nominal_values.txt','w') as f1:
    f1.write('system'+'\t'+
             'mass'+'\t'+
             'age'+'\t'+
             'feh'+'\n')

    for i in range(10):
        system=repr(i+1)

        mass=[]
        age=[]
        feh=[]

        with open('mass_age_teff_sample_'+system+'.txt','r') as f:
            next(f)
            for lines in f:
                x=lines.split()
                try:
                    mass.append(float(x[1]))
                    age.append(float(x[2]))
                    feh.append(float(x[3]))

                except:continue

        mass_mean=corner.quantile(mass,[0.5])[0]
        age_mean=corner.quantile(age,[0.5])[0]
        feh_mean=corner.quantile(feh,[0.5])[0]

        f1.write(system+'\t'+
                 repr(mass_mean)+'\t'+
                 repr(age_mean)+'\t'+
                 repr(feh_mean)+'\n')
