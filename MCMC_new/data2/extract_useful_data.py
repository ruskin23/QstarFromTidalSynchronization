#Extract KIC and data with period less than 10, temperatures more than 5500 and
#less than 6200

def create_array(lines):
    array=[]
    xd = lines.split()
    for value in xd:
        array.append(value)
    return array

#fname = 'kepler_2167_B_data.txt'
fname ='kepler_2167_B_data.txt'

KIC=[]
f1 = open('KIC_for_simbad.txt','w')

with open(fname,'r') as f:
    for line in f:
        if line[0]=='#':continue
        xd = line.split()

        xi=xd[0]
        l=len(xi)
        if xi[0]=='0':a=1
        else:a=0
        y=''.join(xi[k] for k in range(a,l-3))
        xd[0]='KIC'+y

        if float(xd[2])<10.0:
            print(xd[2])
            f1.write(xd[0]+'\n')



f1.close()
