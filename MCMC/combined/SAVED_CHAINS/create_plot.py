import matplotlib.pyplot as plt
import numpy

good_systems=['109','8','123','106','47','84','12','28','70']

mean_values=[]
tidal_days=[]
with open('limits.txt','r') as f:
    next(f)

    for lines in f:
        x=lines.split()
        values=[float(x[1]),float(x[2]),float(x[3])]
        if float(x[6])==0:
            tf=[float(x[7]),float(x[7]),float(x[7])]
            c='b'
            plt.semilogx(tf,values,'-o',color=c)
            x_direct=1
            y_direct=0
            plt.quiver(float(x[7]),float(x[1]),x_direct,y_direct,color='k',width=0.002)

            continue
        if x[0] in good_systems:c='r'
        else:c='g'
        tf=[2*numpy.pi/float(x[6]),2*numpy.pi/float(x[6]),2*numpy.pi/float(x[6])]
        plt.semilogx(tf,values,'-o',color=c)

plt.show()
