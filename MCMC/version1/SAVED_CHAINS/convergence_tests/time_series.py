import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from ks_tests import _fill_parameters

def continuity_check(it,filename):
    for i,j in zip(it,it[1:]):
        if j<i:
            print(filename)
            print(i)
    #print(all(i<j for i,j in zip(it,it[1:])))

#pp = PdfPages('chains_mcmc.pdf')

# plt.rc('xtick',labelsize=30)
# plt.rc('ytick',labelsize=30)

#s=['8','43','36','109','70','47','86','88','93','123','95','106','79','84','25','12','50','28','13']
s=['85', '73', '76', '96', '92', '81', '80', '36', '93', '83', '84', '94', '32', '79', '106', '123', '50', '47', '39', '56', '126', '54', '109', '44', '48', '17', '70', '8', '12', '88', '67', '20', '95', '25', '57', '137', '120', '86', '43', '28', '13']
#s=['1']
clusters=['ganymede','stampede']

parameters=['iteration','porb','eccentricity','wdisk','logQ','mass','age','feh']

for system in s:
    print('\nSystem = ',system)
    # fig,axs=plt.subplots(4,3,figsize=(50,50))
    CHAINS=[]
    ITERATIONS=[]

    for c in clusters:
        for i in range(5):
            n=str(i+1)
            it=[]
            values=[]
            filename='../'+c+'/MCMC_'+system+'/rejected_parameters_'+n+'.txt'
            filled_file=_fill_parameters(filename)
            with open(filled_file,'r') as f:
                #next(f)
                for lines in f:
                    x=lines.split()
                    it.append(int(x[0]))
                    values.append(float(x[4]))
                    #values.append(float(x[parameters.index('logQ')]))
            continuity_check(it,filename)
            ITERATIONS.append(it)
            CHAINS.append(values)
        
#     for i in range(4):
#         for j in range(3):
            
#             if i==3 and j>0:axs[i,j].axis('off')
#             else:axs[i,j].plot(ITERATIONS[3*i+j],CHAINS[3*i+j])
#     plt.suptitle('System '+ system, fontsize=50)
#     plt.savefig(pp,format='pdf')

# pp.close()