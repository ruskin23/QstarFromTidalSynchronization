fname1 ='Kj_2017_ecc_teff_data.txt'
fname2 = 'Kj_2017_period_data.txt'
fout = 'Kj_2017_combined_data.txt'

f3 = open(fout,'w')
f3.write("KIC" + '\t' + 'Teff1' + '\t' + 'Teff2' + '\t' + 'q' + '\t' + 'e' + '\t' + 'Per' + '\n')

with open(fname1,'r') as f1, open(fname2,'r') as f2:
    for line1,line2 in zip(f1,f2):
        d1=line1.split()
        d2=line2.split()
        if 'B' not in d1:
            if float(d2[1]) < 12.0 and float(d1[1])>4800 and float(d1[1])<6200:
                f3.write(d1[0] + '\t' + d1[1] + '\t' + d1[2] + '\t' + d1[3] + '\t' + d1[6] + '\t' + d2[1] + '\n')



f3.close()

