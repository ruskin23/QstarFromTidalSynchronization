import re
fname1='Kj_2017_ecc_teff_data.txt'
fname2='Janes_2017_rot_per.txt'

f3=open('cross_1.txt','w')
with open(fname1,'r') as f1:
    for line in f1:
        name = line.split()
        kic=name[0]
        with open(fname2,'r') as f2:
            for line2 in f2:
                if re.search(kic,line2):
                    #print(line2)
                    f3.write(line2)
                    break
                else:continue
            #    for line2 in f2:
        #        if re.search(kic,line):
        #            f3.write(line2)

f3.close()
