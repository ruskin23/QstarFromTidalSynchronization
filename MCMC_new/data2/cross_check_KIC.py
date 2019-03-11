fname1='eccentricity_teff_data.txt'
fname2='kepler_2167_B_data.txt'

def get_kic(x):
    xi=x[0]
    l=len(xi)
    if xi[0]=='0':a=1
    else:a=0
    y=''.join(xi[k] for k in range(a,l-3))
    return y

with open(fname2,'r') as f2:
    for line in f2:
        x=line.split()
        kic2=get_kic(x)
        if x[0][0]!='#':
            if float(x[2])<10.0 and x[3]=='D' and float(x[6])>5200 and float(x[6])<6000:
                with open(fname1,'r') as f1:
                    for line in f1:
                        y=line.split()
                        kic1=y[0]
                        if kic1==kic2 and float(y[2])>5200 and float(y[2])<6000:
                            print(kic1)
                            break
