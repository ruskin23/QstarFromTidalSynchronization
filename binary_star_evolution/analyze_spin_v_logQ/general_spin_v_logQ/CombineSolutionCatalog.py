


with open('Combined.txt','w') as f:
    f.write('system'+'\t'+
            'logQ'+'\t'+
            'pspin'+'\t'+
            'pspin_error'+'\t'+
            'prob'+'\t'+
            'prob_error'+'\t'+
            'eccentricity'+'\t'+
            'eccentricity_error'+'\t'+
            'tidal_period'+'\t'+
            'tidal_frequency'+'\t'+
            'logtidal_period'+'\t'+
            'logtidal_frequency'+'\t'+
            'primary_mass'+'\t'+
            'age'+'\t'+
            'synchronized'+'\n')
    for system in [1,5,8,12,13,18,20,25,27,28,31,32,33,36,37,39,43,44,47,48,49,50,54,57,59,62,67,68,69,70,73,76,79,80,81,83,84,85,86,88,92,93,94,95,96,99,106,109,111,120,123,125,126,128,129,132,137]:
        data=[]
        print('At System = ', system)
        with open('SolutionFileBreaks0.0.txt','r') as fs:
            next(fs)
            for lines_fs in fs:
                x=lines_fs.split()
                if x[0]==str(system):
                    print('Found in Sol file = ',x[0])
                    break
        with open('SpinlogQCatalog_el0.4.txt','r') as fc:
            next(fc)
            for lines_fc in fc:
                y=lines_fc.split()
                if y[0]==str(system):
                    print('Found in Sol file = ',y[0])
                    break
        data.append(x[0])
        data.append(x[1])
        data.append(x[2])
        data.append(y[13])
        data.append(x[3])
        data.append(y[7])
        data.append(y[8])
        data.append(y[9])
        data.append(x[4])
        data.append(x[5])
        data.append(x[6])
        data.append(x[7])
        data.append(x[8])
        data.append(x[9])
        data.append(x[10])
        values='\t'.join(data)
        f.write(values+'\n')


