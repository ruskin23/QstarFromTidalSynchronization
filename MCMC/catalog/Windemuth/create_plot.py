import matplotlib.pyplot as plt
import numpy
import sys

#n=int(sys.argv[1])

#with open('windemuth_stellar_raw.txt','r') as f:
#    next(f)
#
#    for lines in f:
#        x=lines.split()
#
#        nc_KIC=x[0]
#        xvalues=[]
#        with open('catalog_KIC.txt','r') as g:
#            next(g)
#            for glines in g:
#                y=glines.split()
#                cat_KIC=y[1]
#                if cat_KIC==nc_KIC:
#                    Porb=float(y[6])
#                    xvalues=[Porb,Porb,Porb]
#                    break
#
#        if len(xvalues)==0:continue
#
#        if n==4:
#            age=(10**float(x[n]))/1e9
#            age_error_up=(10**float(x[n+1]))/1e9
#            age_error_down=(10**float(x[n+2]))/1e9
#            age_up=age+age_error_up
#            age_down=age-age_error_down
#            parameter=[age,age_up,age_down]
#        else:parameter=[float(x[n]),float(x[n])+float(x[n+1]),float(x[n])+float(x[n+2])]
#        plt.semilogy(xvalues,parameter,color='k')
#        plt.plot(xvalues[0],parameter[0],'x',color='r')


Porb=[]
e=[]
i=0
with open('windemuth_orbital_raw.txt','r') as f:
    next(f)
    next(f)
    for lines in f:
        x=lines.split()
        nc_KIC=x[0]
        #with open('catalog_KIC.txt','r') as g:
        #    next(g)
        #    for glines in g:
        #        y=glines.split()
        #        cat_KIC=y[1]
        #        if cat_KIC==nc_KIC:
        if float(x[1])<30 and (abs(numpy.sqrt((float(x[7])**2) +
                                              (float(x[10])**2))))<0.5:
            Porb.append(float(x[1]))
            e.append(abs(numpy.sqrt((float(x[7])**2) + (float(x[10])**2))))
        #            break
            i=i+1

print(i)
plt.scatter(Porb,e)
plt.xscale('log')
plt.xlabel('Orbital Period (days)')
plt.ylabel('Eccentricity')
plt.savefig('windemuth_cat.pdf')
