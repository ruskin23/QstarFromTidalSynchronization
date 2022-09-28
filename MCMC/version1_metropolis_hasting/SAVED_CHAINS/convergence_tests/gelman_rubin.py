#system 57 has one chain with 3 accepted. so need to do something bout it

import utils
import numpy
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

def reduce_chain(CHAINS,l_chain):
    #make all chains equal length 

    print(f'min len = {min(l_chain)} max len = {max(l_chain)}')
    min_len=min(l_chain)

    CHAINS_R=[]

    for c in CHAINS:
        chain_sliced=c[:min_len]
        # chain_thinned=chain_sliced[::10]
        CHAINS_R.append(chain_sliced)

    # for _,c in enumerate(CHAINS):
    #     r=int(len(c)/len(CHAINS[min_chain_idx]))
    #     c=c[::r]
    #     diff=len(c)-min_len
    #     if diff==0:
    #         # c=c[::20]
    #         CHAINS_R.append(c)
    #         continue
    #     r=int(len(c)/diff)
    #     ridx=[]
    #     for i in range(len(c)):
    #         if i%(r-1)==0:ridx.append(i)
    #         if len(ridx)==diff:break
    #     c=numpy.delete(c,ridx)
    #     # c=c[::20]
    #     CHAINS_R.append(c)
    
    return CHAINS_R

def adjust_chain(system,chain):

    delete_q=[]
    with open('../deletion_checks/final_results.txt','r') as f:
        for lines in f:
            x=lines.split()
            s=x[0].split('_')[0]
            if s==system:
                if x[6]=='True':
                    delete_q.append(float(x[2]))
                else:delete_q.append(8.353237442958045)

    delete_idx=numpy.zeros(1)
    for q in delete_q:
        idx=numpy.where(chain==q)
        delete_idx=numpy.concatenate((delete_idx,idx[0]),axis=None)
    delete_idx=delete_idx.astype(int)
    chain=numpy.delete(chain,delete_idx)
    return chain



def get_chains(system):
    #returns 10 chains of equal length
    CHAINS=[]
    l_chain=[]

    for c in ['ganymede','stampede']:
        for i in range(1,6):
            chain_filename=utils._get_filename(system,c,i)
            filled_filename=utils._fill_parameters(chain_filename)
            chain=utils._get_chain('logQ',filled_filename)
            chain_fixed=adjust_chain(system,chain)
            l_chain.append(len(chain_fixed))
            CHAINS.append(chain_fixed)
            
    C=reduce_chain(CHAINS,l_chain)
    print(len(C[0]))
    # for i in range(10):
    #     CHAINS[i]=CHAINS[i][-min_len:-1]

    return C

if __name__=='__main__':

    with open('gr_statistics.txt','w') as f:
         f.write('System'+'\t'+'N'+'\t'+'B'+'\t'+'W'+'\t'+'R'+'\n')
    #potential_systems=['8','43','36','109','70','47','86','88','93','123','95','106','79','84','25','12','50','28','13']
    # potential_systems=['70', '36' ,'43' ,'88' ,'86', '8' ,'109', '93', '95' ,'47', '1' ,'123', '106' ,'13', '79', '48' ,'84', '120', '25' ,'12' ,'96' ,'50' ,'28', '94' ,'67']
    non_systems=['76','39','54','80','81','92','126']
    s=['1', '8', '12', '13', '17', '20', '25', '28', '32', '36', '39', '43', '44', '47', '48', '50', '54', '56', '67', '70', '73', '76', '79', '80', '81', '83', '84', '85', '86', '88', '92', '93', '94', '95', '96', '106', '109', '120', '123', '126', '137']
    N=[]
    R=[]
    colors=[]
    for i,system in enumerate(s):
        #print(f'\nSystem = {system}')
        chains=get_chains(system)
        n,b,w,r=gr_test(chains)
        N.append(n)
        R.append(r)
        with open('gr_statistics.txt','a') as f:
            f.write(system+'\t'+
                    repr(n)+'\t'+
                    repr(b)+'\t'+
                    repr(w)+'\t'+
                    repr(r)+'\n')
        print(f'System={system} R={r}')
        if system in non_systems:colors.append('green')
        else:colors.append('red')
    with open('gr_statistics.txt','r') as f:
        next(f)
        for lines in f:
            x=lines.split()
            system=x[0]
            N.append(float(x[1]))
            R.append(float(x[4]))
            if system in non_systems:colors.append('green')
            else:colors.append('red')
    plt.scatter(N,R,color=colors)
    plt.hlines(1,min(N),max(N),linestyles='dashed',color='k',label=1.0)
    plt.hlines(1.1,min(N),max(N),linestyles='dashed',color='b',label=1.1)
    plt.hlines(1.2,min(N),max(N),linestyles='dashed',color='m',label=1.2)
    plt.legend()
    plt.title('Gelman-Rubin Test')
    plt.xlabel('Length of each chain')
    plt.ylabel('R')
    plt.show()
    # plt.savefig('gl_min_start.png')
