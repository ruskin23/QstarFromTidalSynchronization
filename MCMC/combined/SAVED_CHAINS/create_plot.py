import matplotlib.pyplot as plt
import numpy

good_systems=['109','8','123','106','47','84','12','28','70']

mean_values=[]
tidal_days=[]
with open('limits.txt','r') as f:
    next(f)
    k=0
    m=0
    for lines in f:
        x=lines.split()
        values=[float(x[1]),float(x[2]),float(x[3])]
        if float(x[6])==0:
            tf=[float(x[8]),float(x[8]),float(x[8])]
            c='b'
            if k==0:
                plt.semilogx(tf,
                         values,
                         '-o',
                         color=c,
                         label='synchronized systems'
                         )
            else:
                plt.semilogx(tf,
                         values,
                         '-o',
                         color=c
                         )
            x_direct=1
            y_direct=0
            plt.quiver(float(x[7]),
                       float(x[1]),
                       x_direct,
                       y_direct,
                       color='k',
                       width=0.002)

            k=k+1
            continue

        else:
            c='r'
            tf=[2*numpy.pi/float(x[6]),2*numpy.pi/float(x[6]),2*numpy.pi/float(x[6])]
            if m==0:
                plt.semilogx(tf,
                             values,
                             '-o',
                             color=c,
                             label='non-synchronized systems')
            else:
                plt.semilogx(tf,
                             values,
                             '-o',
                             color=c)

            m=m+1

plt.xlabel('Tidal Period (days)')
plt.ylabel('logQ')
plt.legend(prop={'size':6})
#plt.show()
plt.savefig('Results.eps')
