import os
directory=os.getcwd()
print(directory)
numbers=[]
for f in os.listdir(directory):
    if f.endswith('.txt'):
        try:
            x=f.split('_')
            if x[0]=='MassAgeFehSamples':
                y=x[1].split('.')
                numbers.append(y[0])
        except:
            continue


print(numbers)
