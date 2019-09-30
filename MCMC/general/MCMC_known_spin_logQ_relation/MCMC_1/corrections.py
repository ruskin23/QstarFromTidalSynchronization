
D=[]
with open('coorected_combined_accepted.txt','w') as f2:
    f2.write('Iteration_step'+'\t'+
             'primary_mass'+'\t'+
             'age'+'\t'+
             'feh'+'\t'+
             'Porb'+'\t'+
             'eccentricity'+'\t'+
             'Wdisk'+'\t'+
             'logQ'+'\t'+
             'Pspin'+'\n')
    with open('combined_accepted.txt','r') as f:
        for lines in f:
            x=lines.split()
            l=len(x)
            f2.write(x[0]+'\t'+
                 x[1]+'\t'+
                 x[3]+'\t'+
                 x[5]+'\t'+
                 x[7]+'\t'+
                 x[9]+'\t'+
                 x[11]+'\t'+
                 x[13]+'\t'+
                 x[19]+'\n')

