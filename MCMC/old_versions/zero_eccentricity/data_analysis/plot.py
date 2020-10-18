import matplotlib.pyplot as plt

itr=[]
logQ=[]


with open('combined_accepted_parameters.txt','r') as f:
    for lines in f:
        x=lines.split()
        itr.append(int(x[0]))
        logQ.append(float(x[6]))

plt.scatter(itr,logQ)
plt.show()
