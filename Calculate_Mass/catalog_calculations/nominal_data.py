

import re


def write_para(para,m_sol,i):
    filename='catalog_'+repr(i)+'_p.txt'
    with open(filename,'w') as f:
        f.write(kic + '\t' + para[1] + '\t' + para[3] + '\t' + para [9] + '\t' + para[7] + '\t' + para[5] + '\t' + para[11] + '\t' + para[13] + '\n')
    if m_sol==2:
        filename='catalog_'+repr(i)+'_s.txt'
    with open(filename,'w') as f:
        f.write(kic + '\t' + para[1] + '\t' + para[3] + '\t' + para [9] + '\t' + para[7] + '\t' + para[5] + '\t' + para[11] + '\t' + para[13] + '\n')

i=1
with open('mass_upto_0.2.txt','r') as f1:
    next(f1)
    for line1 in f1:
        data=line1.split()
        kic=data[0]
        print(kic)
        with open('KIC_catalog.txt','r') as f2:
            next(f2)
            for line2 in f2:
                print('checking')
                if re.search(kic,line2):
                    print('found')
                    para=line2.split()
                    if len(data)==7:m_sol=2
                    else:m_sol=1
                    write_para(para,m_sol,i)
                    i=i+1
                    break
                else:
                    continue
