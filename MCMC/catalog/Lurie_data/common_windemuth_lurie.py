

Lurie=[]
with open('binaries_with_p1.txt','r') as f:
    for lines in f:
        x=lines.split()
        Lurie.append(x[0])

n=0
with open('Windemuth_combined_filtered.txt','r') as f:
    for lines in f:
        x=lines.split()
        if x[1] in Lurie:
            n=n+1
print(n)
