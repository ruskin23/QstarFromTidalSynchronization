import os
import numpy

#systems=[ '18', '23', '27', '33', '35', '37', '49', '59', '62']
systems=['2', '5', '14', '18', '23', '27' ,'33' ,'35', '37', '49', '59', '62', '63', '65', '68', '69', '77', '78', '82', '99', '100', '111', '117', '119', '125', '128', '129', '132', '133', '134' ,'138', '142']
for system in systems:

    print(system)
    directory='System_'+system

    dataset=numpy.zeros([7,11])
    dataset[1:7,0]=[k for k in range(5,11)]
    percentile_values=[1,2,3,4,5,10,20,30,40,50]
    dataset[0,1:12]=[k for k in percentile_values]

    catalog_file='../../SpinlogQCatalog_el0.4.txt'
    with open(catalog_file,'r')as f:
        next(f)
        for lines in f:
            x=lines.split()
            at_system=x[0]
            if at_system==system:
                Spin_value=float(x[12])
                break
    dataset[0,0]=system

    for filename in os.listdir(directory):
        print(filename)
        percentile = int((filename.split('_')[1]).split('.')[0])
        if percentile//10==0:percentile_index=percentile
        else:percentile_index=(percentile//10)+5
        with open(os.path.join(directory,filename),'r') as f:
            next(f)
            for i,lines in enumerate(f):
                x=lines.split()
                dataset[i+1,percentile_index]=float(x[1])-Spin_value


    datafile='data_file_'+system+'.txt'
    with open(datafile,'w') as f:
        for row in dataset:
            v=[repr(r) for r in row]
            f.write('\t'.join(v)+'\n')
        f.write('\n')
