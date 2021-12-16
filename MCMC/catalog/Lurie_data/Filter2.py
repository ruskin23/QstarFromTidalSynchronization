#filter to select binaries with Porb<60

fo=open('filter2.txt','w')
with open('binaries_with_p1.txt','r') as f:
    for lines in f:
        x=lines.split()
        if float(x[1])<60:
            fo.write(lines)
fo.close()
