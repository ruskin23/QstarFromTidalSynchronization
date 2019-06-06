import re

fname1='cross_2.txt'
fname2='Kj_2017_ecc_teff_data.txt'
fname3='Mathur_2017_data.txt'
fname4='KEBC_period_data.txt'

with open('catalog_KIC.txt','w') as f:
    f.write('KIC' + '\t' + 'teff' + '\t' + 'teff_e' + '\t' + 'feh' + '\t' +
            'feh_e' + '\t' + 'Porb' + '\t' + 'Porb_e' + '\t' + 'e' + '\t' +
            'e_e' + '\t' + 'logg' + '\t' + 'logg_e'  + '\t' + 'Pspin' + '\t' +
            'Pspin_e' + '\t' + 'q' + '\n')


KIC=[]

teff=[]
teff_e=[]

feh=[]
feh_e=[]

Porb=[]
Porb_e=[]

e=[]
e_e=[]

logg=[]
logg_e=[]

Pspin=[]
Pspin_e=[]

q=[]

fc=open('catalog_KIC.txt','a')
with open(fname1,'r') as f1: #Lurie:
    for line in f1:
        values_lurie = line.split()
        kic=values_lurie[0]
        porb=float(values_lurie[1])
        pspin=float(values_lurie[3])
        pspin_max=float(values_lurie[4])
        with open(fname2,'r') as f2:#Kj
            for line2 in f2:
                if re.search(kic,line2):
                    with open(fname3,'r') as f3:#MAth
                        for line3 in f3:
                            if re.search(kic,line3):
                                values_mathur=line3.split()
                                values_kj=line2.split()

                                KIC.append(int(kic))
                                teff.append(float(values_mathur[1]))
                                teff_e.append(float(values_mathur[2]))

                                feh.append(float(values_mathur[7]))
                                feh_e.append(float(values_mathur[8]))

                                logg.append(float(values_mathur[4]))
                                logg_e.append(float(values_mathur[5]))

                                Porb.append(porb)

                                evalue=float(values_kj[6])
                                e.append(evalue)
                                if evalue>0.1:e_e.append(0.001)
                                if evalue<0.1 and evalue>0.01:e_e.append(0.01)
                                if evalue<0.01:e_e.append(0.1)

                                Pspin.append(pspin)
                                spin_error=pspin_max-pspin
                                Pspin_e.append(spin_error)

                                if values_kj[2]=='B':q.append(float(values_kj[4]))
                                else:q.append(float(values_kj[3]))


        with open(fname4,'r') as f4:
            for line4 in f4:
                if re.search(kic,line4):
                    values=line4.split()
                    Porb_e.append(values[2])


for i in range(len(KIC)-1):
    if Porb[i]<60 and e[i]<0.01:
        fc.write(repr(KIC[i]) + '\t' + repr(teff[i]) + '\t' + repr(teff_e[i]) + '\t' + repr(feh[i]) + '\t' + repr(feh_e[i]) + '\t' + repr(Porb[i]) + '\t' + Porb_e[i] + '\t' + repr(e[i]) + '\t' + repr(e_e[i]) + '\t' + repr(logg[i]) + '\t' + repr(logg_e[i]) + '\t' + repr(Pspin[i]) + '\t' + repr(Pspin_e[i]) + '\t' + repr(q[i]) + '\n')





fc.close()

