import numpy


with open('catalog_l0.4.txt','r') as f:
    for i,lines in enumerate(f):
        x=lines.split()
        print('At System = ',x[0])
        if i==0:
            with open('SpinlogQCatalog_el0.4.txt','w') as fH:
                fH.write('\t'.join(x)+'\n')
                continue
        primary_mass=float(x[15])
        mass_ratio=float(x[14])
        secondary_mass=primary_mass*mass_ratio
        print('PrimaryMass = {}, MassRatio = {}, SecondaryMass = {}'.format(primary_mass,mass_ratio,secondary_mass))
        if numpy.logical_and(secondary_mass>0.4,
                            secondary_mass<1.2):
            with open('SpinlogQCatalog_el0.4.txt','a') as fN:
                values='\t'.join(x)
                fN.write(values+'\n')
