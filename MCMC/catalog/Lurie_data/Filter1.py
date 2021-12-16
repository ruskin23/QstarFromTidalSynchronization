#filter to see how many binaries have P1min

fo=open('binaries_with_p1.txt','w')
with open('sp_binaries.txt','r') as f:
    for lines in f:
        x=lines.split()
        if len(x)>=6 and x[-1]!='b':
            fo.write(lines)
fo.close()
