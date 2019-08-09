
terms_age = ['Age:','Age:|']
age=[]

terms_condition = ['Condition[0]:','Condition[0]:|']
condition=[]

with open('test_data.txt','r') as f:
    for line in f:
        d=line.split()
        if d[0] in terms_age:
            del d[0]
            for item in d:
                if '|' in item:item=item.replace('|','')
                age.append(float(item))
        if len(d)>1:
            print(d)
            a=d[0]+d[1]
            if a in terms_condition:
                for i in range(4):
                    del d[0]
                for item in d:
                    if '|' in item:item=item.replace('|','')
                    condition.append(float(item))


print(age)
print(condition)
