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
de_upper=[]
de_lower=[]
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
        e1=float(x[7])
        e2=float(x[10])
        e_value=abs(numpy.sqrt(e1**2 + e2**2))
        delta_e1_upper=float(x[8])
        delta_e1_lower=float(x[9])
        delta_e2_upper=float(x[11])
        delta_e2_lower=float(x[12])
        e1_r=e1/e_value
        e2_r=e2/e_value
        if float(x[1])<30 and e_value<0.5:
            Porb.append(float(x[1]))
            e.append(e_value)
            de_upper.append((e1_r*delta_e1_upper) + (e2_r*delta_e2_upper))
            de_lower.append((e1_r*delta_e1_lower) + (e2_r*delta_e2_lower))
            i=i+1

print(i)
e=numpy.array(e)
de_lower=numpy.array(de_lower)
de_upper=numpy.array(de_upper)
errors=numpy.vstack([de_lower,de_upper])

plt.scatter(Porb,e)
plt.xscale('log')
plt.xlabel('Orbital Period (days)')
plt.ylabel('Eccentricity')
plt.show()
#plt.savefig('windemuth_cat.pdf')
