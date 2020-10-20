
with open('parameters.txt','w') as f:
    f.write('orbital_period'+'\t'+
            'feh'+'\t'+
            'eccentricity'+'\t'+
            'Wdisk'+'\t'+
            'logQ'+'\t'+
            'primary_mass'+'\t'+
            'age'+'\n')

values=[]
with open('temp.out','r') as f:
    for i,lines in enumerate(f):
        x=lines.split()
        if x[0]=='Parameter':continue
        elif x[0]=='--':
            values.append('\n')
            with open('parameters.txt','a') as g:
                g.write('\t'.join(values))
            values=[]
        else:values.append(x[2])


