import matplotlib.pyplot as plt
import numpy
s=['1', '8', '12', '13', '17', '20', '25', '28', '32', '36', '39', '43', '44', '47', '48', '50', '54', '56', '67', '70', '73', '76', '79', '80', '81', '83', '84', '85', '86', '88', '92', '93', '94', '95', '96', '106', '109', '120', '123', '126', '137']


porb=[]
pspin=[]
r=[]
e=[]
with open('catalog_KIC.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        #if x[0] in s:
            #if float(x[8])<0.1:
            #    print(f'{x[0]} {float(x[6])} {float(x[8])}')
        porb.append(float(x[6]))
        pspin.append(float(x[11]))
        #if float(x[5])/float(x[11])<6:
        r.append(float(x[5])/float(x[11]))
        e.append(float(x[8]))
        #if float(x[5])/float(x[11])>1.2 and float(x[7])>0.1 and float(x[7])<0.3:
        #print(lines)


e=numpy.array(e)
porb=numpy.array(porb)
pspin=numpy.array(pspin)
r=numpy.array(r)

plt.scatter(porb,e)
plt.xscale('log')
plt.show()


