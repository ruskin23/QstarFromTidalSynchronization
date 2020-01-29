

n=[]
with open('NewAgeSystems.txt','r') as f:
    for lines in f:
        x=lines.split()
        n.append(str(x[2]))
print(n)
print(' '.join(n))
