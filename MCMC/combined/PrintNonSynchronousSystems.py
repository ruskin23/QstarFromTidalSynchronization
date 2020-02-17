

s=[]
with open('SolutionFileBreaks0.0.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        if x[10]=='None':
            s.append(x[0])


print(' '.join(s))
