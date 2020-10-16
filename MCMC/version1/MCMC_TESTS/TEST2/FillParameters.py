



with open('AccetedParameters.txt','w') as f1:
    with open('accepted_parameters_1.txt','r') as f2:
        for i,lines in enumerate(f2):
            if i==0:
                f1.write(lines)
                continue
            if i>0:

                print('\nAt i = ', i)
                x=lines.split()
                current_state=x[1:-1]
                IterationNumber=float(x[0])
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

