import sys


system=sys.argv[1]
instance=sys.argv[2]
system_filename='MCMC_'+system+'/accepted_parameters_'+instance+'.txt'
with open('AcceptedParameters.txt','a') as f1:
    with open(system_filename,'r') as f2:
        for i,lines in enumerate(f2):
            if i==0:
                continue
            if i>0:

                print('\nAt i = ', i)
                x=lines.split()
                current_state=x[1:-1]
                IterationNumber=int(x[0])
                if i>1:
                    if previous_state[0]==x[0]:continue
                if IterationNumber==1:
                    print('First Iteration')
                    Step=1
                    previous_state=current_state
                    f1.write(lines)
                    continue
                else:
                    Step=Step+1
                    if IterationNumber!=Step:
                        while True:
                            if Step!=IterationNumber:
                                print('Step = ',Step)
                                print('IterationNumber = ',IterationNumber)
                                state=[str(Step)]+previous_state
                                Parameters='\t'.join(state)
                                f1.write(Parameters+'\n')
                                Step=Step+1
                            else:
                                state=[str(Step)]+current_state
                                Parameters='\t'.join(state)
                                f1.write(Parameters+'\n')
                                previous_state=current_state
                                break
                    else:
                        state=[str(Step)] + current_state
                        Parameters='\t'.join(state)
                        f1.write(Parameters+'\n')
                        previous_state=current_state

