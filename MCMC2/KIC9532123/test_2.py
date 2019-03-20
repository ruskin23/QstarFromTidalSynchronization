import numpy
def fun(x):
    if x>5: raise AssertionError
    y=42
    print('am i stupid or this thing does not work = ', y)
    print('value of x in fun = ',x)

    return x
x=126
while True:
    try :
        print(fun(x))
        break
    except AssertionError:
        print('value ', x)
        x = x-1
        continue

