import matplotlib.pyplot as plt
import numpy

good_systems=['70','43','36','88','86','93','8','95','109','1','47','123','48','13','106','96','79','94','67','84','120','25','83','12','50','57','28','17','44','85']


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
plt.show()
#plt.savefig('Results.eps')

