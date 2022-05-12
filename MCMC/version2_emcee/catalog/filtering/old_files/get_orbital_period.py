fo=open('orbital_periods.txt','w')
with open('windemuth_orbital_raw.txt','r') as f:
    next(f)
    next(f)
    for lines in f:
        x=lines.split()
        fo.write(x[0]+'\t'+x[1]+'\n')