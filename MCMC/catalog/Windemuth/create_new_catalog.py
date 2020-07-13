import numpy


with open('WindeCatalog.txt','w') as f:
    f.write('system'+'\t'+
            'KIC'+'\t'+
            'Porb'+'\t'+
            'Porb_sigma'+'\t'+
            'ecosw'+'\t'+
            'ecosw_sigma'+'\t'+
            'esinw'+'\t'+
            'esimw_sigma'+'\t'+
            'M1'+'\t'+
            'M1_sigma'+'\t'+
            'M2'+'\t'+
            'M2_sigma'+'\t'+
            'age(log10 yr)'+'\t'+
            'age_sigma'+'\t'+
            'feh'+'\t'+
            'feh_sigma'+'\t'+
            'Pspin'+'\t'+
            'Pspin_sigma'+'\n')

with open('catalog_KIC.txt','r') as f:
    next(f)
    for flines in f:
        x=flines.split()
        cat_KIC=x[1]

        with open('windemuth_stellar_raw.txt','r') as g:
            for glines in g:
                y=glines.split()
                w_KIC=y[0]
                if cat_KIC==w_KIC:

                    with open('windemuth_orbital_raw.txt','r') as h:
                        for hlines in h:
                            z=hlines.split()
                            wo_KIC=z[0]
                            if w_KIC==wo_KIC:



                                viable=True
                                #if float(z[1])>40:viable=False
                                #if numpy.sqrt((float(z[7])**2) + (float(z[10])**2))>0.50:viable=False
                                if float(y[7])>1.2 or float(y[7])<0.4:viable=False
                                if float(y[10])>1.2 or float(y[10])<0.4:viable=False
                                if float(y[1])<-1.014 or float(y[1])>0.537:viable=False

                                if viable==True:
                                    ecc=numpy.sqrt((float(z[7])**2) + (float(z[10])**2))
                                    ecc_error=((float(z[7])*float(z[8])) +
                                               (float(z[10])*float(z[11])))/ecc
                                    with open('WindeCatalog.txt','a') as n:
                                        n.write(x[0]+'\t'+
                                                x[1]+'\t'+
                                                z[1]+'\t'+
                                                repr(abs(float(z[2])))+'\t'+
                                                repr(float(z[7]))+'\t'+
                                                repr(abs(float(z[8])))+'\t'+
                                                repr(float(z[10]))+'\t'+
                                                repr(abs(float(z[11])))+'\t'+
                                                y[7]+'\t'+
                                                repr(abs(float(y[8])))+'\t'+
                                                y[10]+'\t'+
                                                repr(abs(float(y[11])))+'\t'+
                                                repr(float(y[4]))+'\t'+
                                                repr(float(y[5]))+'\t'+
                                                y[1]+'\t'+
                                                repr(abs(float(y[2])))+'\t'+
                                                x[12]+'\t'+
                                                x[13]+'\n')
