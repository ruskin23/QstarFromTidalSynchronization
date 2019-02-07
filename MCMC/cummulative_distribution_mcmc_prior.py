import csv
import numpy as np
import scipy
from scipy.stats import norm

import matplotlib.pyplot as plt

iterations = []
value =[]
with open('accepted_test_1.txt', 'r') as f:
    reader = csv.reader(f, dialect='excel-tab')
    for line in reader:
        iterations.append(line[0:1])
        value.append(line[0:4])


size = len(iterations)

for index in range(size):
    for i in iterations[index]:
        l = int(i)
        iterations[index] = l

    for j in value[index]:
        m = float(j)
        value[index] = m



total_iterations = iterations[size - 1]

#check_total_iterations = iterations[10]

check = np.zeros(total_iterations)
value_array = np.zeros(total_iterations)

print(total_iterations)
print(len(value_array))

for i in range(size - 1):
    value_array[iterations[i]] = value[i]

for i,v in  enumerate(value_array):
    if v!=0: 
        k = i 
        break

while True:

    if k > len(value_array) - 1 : break    
    if value_array[k]==0: value_array[k] = value_array[k-1]
    k = k+1


value_array.sort()

for x in value_array:
    print(x)
#yaxis = numpy.linspace(0,value_array[len(value_array) - 1 ],len(value_array) - 1 )


#dist = scipy.stats.norm.cdf(value_array)

#dist_norm = np.random.normal(5922,100)


#mu, sigma = 0, 0.1 # mean and standard deviation
#s = np.random.normal(mu, sigma, 1000)

#plt.plot(s)
#plt.show()

