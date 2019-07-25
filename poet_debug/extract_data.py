#Exctract data from stopping condition lines.
import matplotlib.pyplot as plt

fname ='last_stopping_condition.txt'

with open(fname,'r') as f:
    for i,line in enumerate(f):
        value = line.split()
        if i==0:age=value
        if i==1:condition0=value
        if i==2:condition1=value
        if i==3:condition2=value
        if i==4:condition3=value
        if i==5:condition4=value
        if i==6:condition5=value
        if i==7:condition6=value

remove_con0=['Condition[', '0]:', '\\', 'ze>']
remove_con1=['Condition[', '1]:', '\\', 'ze>']
remove_con2=['Condition[', '2]:', '\\', 'ze>']
remove_con3=['Condition[', '3]:', '/', 'ze>']
remove_con4=['Condition[', '4]:', '/', 'ze>']
remove_con5=['Condition[', '5]:', '\\', 'ze>']
remove_con6=['Condition[', '6]:', '/', 'ze>']
remove_age=['Age:']

for items in remove_age:
    age.remove(items)
for items in remove_con0:
    condition0.remove(items)
for items in remove_con1:
    condition1.remove(items)
for items in remove_con2:
    condition2.remove(items)
for items in remove_con3:
    condition3.remove(items)
for items in remove_con4:
    condition4.remove(items)
for items in remove_con5:
    condition5.remove(items)
for items in remove_con6:
    condition6.remove(items)




plt.plot(age,condition2)
plt.show()
