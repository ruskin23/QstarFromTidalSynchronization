
import matplotlib.pyplot as plt
age=[]
with open('test.txt','r') as f:
    for lines in f:
        x=lines.split()
        age.append(float(x[3]))


plt.plot(age)
plt.show()


