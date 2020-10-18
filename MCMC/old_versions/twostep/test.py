from attributedict.collections import AttributeDict

parameters=dict(a=dict(value=10,sigma=12),
                b=dict(value=13,sigma=15),
                c=123123)


#p=AttributeDict({'foo': {'bar': [1, 2, 3]}})
p=AttributeDict(parameters)
print(p.c)
