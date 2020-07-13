import sys

system = sys.argv[1]
cluster=['ganymede', 'stampede']
instances=['1','2','3','4','5']

print('\nSystem = ',system)

for c in cluster:
    print('cluster = ',c)
    cluster_directory='/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/combined/SAVED_CHAINS/'+c
    system_directory=cluster_directory+'/MCMC_'+system
    for i in instances:
        print('instance = ',i)
        accepted_filename=system_directory+'/accepted_parameters_'+i+'.txt'
        with open(accepted_filename,'r') as f:
            next(f)
            next(f)
            for lines in f:
                x=lines.split()
                with open('extracted.txt','a') as fe:
                    fe.write(x[1]+'\t'+
                             x[2]+'\t'+
                             x[3]+'\t'+
                             x[4]+'\t'+
                             x[5]+'\t'+
                             x[6]+'\t'+
                             x[7]+'\t'+
                             x[8]+'\t'+
                             x[11]+'\t'+
                             x[12]+'\n')
