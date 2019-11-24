import numpy as np
import matplotlib.pyplot as plt
import corner

index=[]
with open('catalog_KIC.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        index.append(x[0])



with open('nominal_values.txt','w') as f1:
    f1.write('system'+'\t'+
             'mass'+'\t'+
             'mass_sigma'+'\t'+
             'age'+'\t'+
             'age_sigma'+'\t'+
             'feh'+'\t'+
             'feh_sigma'+'\n')

    for system in index:

        mass=[]
        age=[]
        feh=[]

        if system in ['111','139']:continue

        with open('MassAgeFehSamples_'+system+'.txt','r') as f:
            next(f)
            for lines in f:
                x=lines.split()
                try:
                    mass.append(float(x[1]))
                    age.append(float(x[2]))
                    feh.append(float(x[3]))

                except:continue
        mass.sort()
        age.sort()
        feh.sort()
        mass_mean=corner.quantile(mass,[0.5])[0]
        age_mean=corner.quantile(age,[0.5])[0]
        feh_mean=corner.quantile(feh,[0.5])[0]

        mass_sigma=abs(corner.quantile(mass,[0.6827])[0]-max(mass))
        age_sigma=abs(corner.quantile(age,[0.6827])[0]-max(age))
        feh_sigma=abs(corner.quantile(feh,[0.6827])[0]-max(feh))

        f1.write(system+'\t'+
                 repr(mass_mean)+'\t'+
                 repr(mass_sigma)+'\t'+
                 repr(age_mean)+'\t'+
                 repr(age_sigma)+'\t'+
                 repr(feh_mean)+'\t'+
                 repr(feh_sigma)+'\n')
