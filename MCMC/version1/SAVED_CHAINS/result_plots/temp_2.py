import matplotlib.pyplot as plt
import numpy

s=['1', '8', '12', '13', '17', '20', '25', '28', '32', '36', '39', '43', '44', '47', '48', '50', '54', '56', '67', '70', '73', '76', '79', '80', '81', '83', '84', '85', '86', '88', '92', '93', '94', '95', '96', '106', '109', '120', '123', '126', '137']

s1=[]
with open('../SpinlogQCatalog_el0.4.txt','r') as f:
    for lines in f:
        x=lines.split()
        if x[0] not in s:
            print(f'{x[0]} {x[1]} ')
            s1.append(x[1])
with open('temp_5.txt','r') as f:
    for lines in f:
        x=lines.split()
        if x[0] in s1:
            print(x[0])

# E_p=[]
# with open('E_p.txt','r') as f:
#     next(f)
#     for lines in f:
#         x=lines.split()
#         E_p.append(float(x[1]))


# E_p.sort()
# E_p=numpy.array(E_p)
# x=numpy.arange(len(E_p))


# fig=plt.figure()
# ax=plt.gca()
# ax.scatter(x,E_p)
# ax.set_yscale('log')
# plt.xlabel('Index')
# plt.ylabel(r'$\log_{10}{E(p)}$')
# plt.savefig('E_p.png')