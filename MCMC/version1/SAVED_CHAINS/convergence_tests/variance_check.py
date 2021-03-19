import numpy
import utils
import pickle
from scipy.stats import levene,bartlett,fligner,chi2
import matplotlib.pyplot as plt

def degree_of_freedoms(N,M,s_array,x_mean_array,total_mean,B,V):

    f1 = (((N-1)/N)**2)*(1/M)*numpy.var(s_array)
    f2 = (((M+1)/(M*N))**2)*(2/(M-1))*B*B
    f3 = 2*((M+1)*(N+1))/(M*N*N)
    # X=numpy.stack((s_array,x_mean_array**2),axis=0)
    f4 = (N/M)*(numpy.cov(s_array,x_mean_array**2)[0][1]-2*total_mean*numpy.cov(s_array,x_mean_array)[0][1])
    print(f'f1 = {f1} f2 = {f2} f3 = {f3} f4 = {f4}')

    var_V = f1 + f2 + f3*f4
    df = (2*V*V)/(var_V)
    return df

def gr_test(chains):
 
    N=len(chains[0])
    M=10
    chains=numpy.array(chains)
    total_mean=numpy.mean(chains)

    #between chain variance
    B=0
    x_mean_array=[]
    for c in chains:
        B=B+((numpy.mean(c)-total_mean)**2)
        x_mean_array.append(numpy.mean(c))
    x_mean_array=numpy.array(x_mean_array)
    B=B*(N/(M-1))
    # B=N*numpy.mean((numpy.mean(chains,axis=1)-total_mean)**2,)*(M/(M-1))

    #within chain variance
    W=0
    s_array=[]
    for c in chains:
        W=W+(numpy.var(c))
        s_array.append(numpy.var(c))
    s_array=numpy.array(s_array)
    W=W/M
    # W=numpy.mean(numpy.var(chains,axis=1))

    #target variance
    # var_target=(((N-1)/N)*W) + (((M+1)/(M*N))*B)
    var_target=((N-1)/N)*W + B/N
    V = var_target + B/(M*N)

    #correction
    # df=degree_of_freedoms(N,M,s_array,x_mean_array,total_mean,B,V)

    #GR statistics
    R=V/W
    R=numpy.sqrt(R)
    # R=numpy.sqrt(R*(df/(df-2)))

    return N,B,W,R


def get_chain(system):

    with open('../complete_chains.pickle','rb') as f:
        D=pickle.load(f)
    for system_name,params in D.items():
        if system_name==system:
            for name,values in params.items():
                if name=='logQ':main_chain=values

    return main_chain

def divide_chain(system,div):

    main_chain=numpy.zeros(1)    
    chains=[[1] for k in range(10)]
    for clust in ['ganymede','stampede']:

        for i in range(1,6):

            chain_filename=utils._get_filename(system,clust,i)
            filled_filename=utils._fill_parameters(chain_filename)
            chain=utils._get_chain('logQ',filled_filename)
            chain_fixed=utils.adjust_chain(system,chain)
            max_len=len(chain_fixed)
            end_len=int(div*max_len)
            chain_fixed=chain_fixed[:end_len]
            chain_split=numpy.array_split(chain_fixed,10)
            for k,c in enumerate(chain_split):
               for v in c:
                   chains[k].append(v)
    #         if chain_filename!=f'ganymede/MCMC_{system}/accepted_parameters_1.txt':
    #             start_point=chain_fixed[0]
    #             count_repeat=numpy.count_nonzero(chain_fixed==start_point)
    #             chain_fixed=chain_fixed[count_repeat:]
    #         main_chain=numpy.concatenate((main_chain,chain_fixed),axis=None)
    
    # main_chain=main_chain[1:]
    
    return chains
    # main_chain=get_chain(system)
    # return numpy.array_split(main_chain,20)
   

def variance_test(splits,test_name):

    if test_name=='levene':
        test_obj=levene
    elif test_name=='bartlett':
        test_obj=bartlett
    elif test_name=='fligner':
        test_obj=fligner

    stat,p=test_obj(splits[0],
                    splits[1],
                    splits[2],
                    splits[3],
                    splits[4],
                    splits[5],
                    splits[6],
                    splits[7],
                    splits[8],
                    splits[9],
                    splits[10],
                    splits[11],
                    splits[12],
                    splits[13],
                    splits[14],
                    splits[15],
                    splits[16],
                    splits[17],
                    splits[18],
                    splits[19]
                    )

    return stat,p

s=['1', '8', '12', '13', '17', '20', '25', '28', '32', '36', '39', '43', '44', '47', '48', '50', '54', '56', '67', '70', '73', '76', '79', '80', '81', '83', '84', '85', '86', '88', '92', '93', '94', '95', '96', '106', '109', '120', '123', '126', '137']
s=['137']
stats=[]
L=[]
# with open('new_gr_test.txt','w',1) as f:
#     f.write('KIC\tR\n')
for system in s:
    for d in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]:
        splits=divide_chain(system,d)
        l_chain=[]
        chains=[]
        for ch in splits:
            l_chain.append(len(ch))
        min_l=min(l_chain)
        for ch in splits:
            chains.append(ch[:min_l])


        N,B,W,R=gr_test(chains)
        stats.append(R)
        L.append(numpy.log10(min_l))
        with open('../SpinlogQCatalog_el0.4.txt','r') as f:
            next(f)
            for lines in f:
                X=lines.split()
                if X[0]==system:
                    KIC=X[1]
                    break
            # f.write(KIC+'\t'+repr(R)+'\n')
        print(f'{KIC} {B} {W} {R}')

# plt.scatter(L,stats)
# plt.axhline(1.1,color='green',label='1.1')
# plt.axhline(1.2,color='red',label='1.2')
# plt.savefig('new_gr_test.png')


#derive constraints only from e>0.1 and the 2 with e~0.001
