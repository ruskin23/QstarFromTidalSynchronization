


with open('AcceptedParameters.txt','r') as f:
    line=list(f)[-1]
    x=line.split()
    print(int(x[0]))
