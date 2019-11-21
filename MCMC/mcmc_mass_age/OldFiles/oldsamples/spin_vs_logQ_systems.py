import numpy

index=[]
with open('catalog_KIC.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        index.append(x[0])


with open('spin_vs_logQ_systems_0.2_0.4.txt','w') as fs:
    fs.write('Head'+'\n')
    for system in index:

        if system in ['111','139']:continue

        with open('catalog_KIC.txt','r') as f:
            for lines in f:
                x=lines.split()
                if x[0]==system:
                    cat_line=lines
                    e=float(x[8])
                    q=float(x[14])
                    break

        with open('Corrected_Nominal.txt') as f:
            for lines in f:
                y=lines.split()
                if y[0]==system:
                    mass_line=lines
                    y=lines.split()
                    primary_mass=float(y[1])
                    break


        secondary_mass=primary_mass*q

        if numpy.logical_and(numpy.logical_and(secondary_mass>0.4,
                             secondary_mass<1.2),
                             e>0.2):
            print('q = ',q)
            print('primary_mass = ', primary_mass)
            print('Secondary Mass = ',secondary_mass)

            z=[y[i] for i in range(1,4)]
            z=x+z
            full_line='\t'.join(z)
            fs.write(full_line+'\n')



