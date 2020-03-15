



with open('AccetedParameters.txt','w') as f1:
    with open('combined_accepted.txt','r') as f2:
        for i,lines in enumerate(f2):
            if i==0:
                f1.write(lines)
                continue
            if i>0:

                print('\nAt i = ', i)
                x=lines.split()
                IterationNumber=float(x[0])

                if IterationNumber==1:
                    print('First Iteration')
                    PreviousIterationStep=1
                    f1.write(lines)
                    continue
                else:

                    if IterationNumber!=PreviousIterationStep:
                        while IterationNumber!=PreviousIterationStep:
                            print('IterationNumber = ',IterationNumber)
                            print('PreviousIterationStep = ', PreviousIterationStep)
                            x[0]=str(PreviousIterationStep+1)
                            Parameters='\t'.join(x)
                            f1.write(Parameters+'\n')
                            PreviousIterationStep=PreviousIterationStep+1
                    else:
                        Parameters='\t'.join(x)
                        f1.write(Parameters+'\n')

