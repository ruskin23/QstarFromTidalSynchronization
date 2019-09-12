data = {'age':[2],
        'pm':[5]}

data['age'].append(7)
data['pm'].append(99)

print(len(data['age']))

for i in range(len(data['age'])):
    print('inst = ',i)
    print(data['age'][i])
    print(data['pm'][i])


print(data)
