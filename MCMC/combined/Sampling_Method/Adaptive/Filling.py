import sys


Step=1
case=sys.argv[1]
i=sys.argv[2]
with open(case+'/AcceptedParameters.txt','w') as f:
    with open(case+'/accepted_parameters_'+i+'.txt','r') as fa:
        next(fa)
        for i,lines in enumerate(fa):
            x=lines.split()
            current_iteration=int(x[0])
            current_parameters=x[1:-1]
            if i==0:
                f.write(lines)
                Step=Step+1
                previous_parameters=current_parameters
                continue
            if i>0:
                if current_iteration!=Step:
                    while current_iteration!=Step:
                        step=[str(Step)]
                        p=step+previous_parameters
                        parameters='\t'.join(p)
                        f.write(parameters+'\n')
                        Step=Step+1
                else:
                    f.write(lines)

            previous_parameters=current_parameters
