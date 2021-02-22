import pickle
import matplotlib.pyplot as plt
import numpy


x=numpy.linspace(5,12,100000)
with open('all_pdf_data.pickle','rb') as f:
    D=pickle.load(f)

s=['85', '73', '76', '96', '92', '81', '80', '36', '93', '83', '84', '94', '32', '79', '106', '123', '50', '47', '39', '56', '126', '54', '109', '44', '48', '17', '70', '8', '12', '88', '67', '20', '95', '25', '137', '120', '86', '43', '28', '13']
s=numpy.reshape(numpy.array(s),(8,5))

M=D['M']/max(D['M'])
y=D['12']/max(D['12'])

plt.plot(x,M)
plt.plot(x,y)
plt.show()

# fig,axs=plt.subplots(8,5,sharex='col',sharey='row',gridspec_kw={'hspace': 0, 'wspace': 0})



# for i in range(8):
#     for j in range(5):
#         print(s[i][j])
#         print(len(x))
#         print(len(D[s[i][j]]))
#         axs[i,j].plot(x,D[s[i][j]])
#         axs[i,j].plot(x,M)
#         axs[i,j].tick_params(axis='x', labelsize= 3)
#         axs[i,j].tick_params(axis='y', labelsize= 3)
#         with open('../../SpinlogQCatalog_el0.4.txt','r') as f:
#             next(f)
#             for lines in f:
#                 y=lines.split()
#                 if y[0]==s[i][j]:
#                     KIC=y[1]
#                     break
#         # axs[i,j].set_title(KIC)

# for ax in axs.flat:
#     ax.label_outer()


# plt.savefig('all_pdf.png')