#system 57 has one chain with 3 accepted. so need to do something bout it

import utils
import numpy
import matplotlib.pyplot as plt


def gr_test(chains):

    N=len(chains[0])
    M=10
    chains=numpy.array(chains)
    total_mean=numpy.mean(chains)

    #between chain variance
    B=N*numpy.mean((numpy.mean(chains,axis=1)-total_mean)**2,)*(M/(M-1))

    #within chain variance
    W=numpy.mean(numpy.var(chains,axis=1))

    #target variance
    var_target=((N-1)/N)*W + B/N

    #GR statistics
    R=numpy.sqrt(var_target/W)

    return N,B,W,R

def reduce_chain(CHAINS,l_chain):
    #make all chains equal length 

    print(f'min len = {min(l_chain)} max len = {max(l_chain)}')
    min_len=min(l_chain)
    max_len=max(l_chain)
    min_chain_idx=l_chain.index(min_len)

    CHAINS_R=[]

    for c in CHAINS:
        chain_sliced=c[:min_len]
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
    plt.show()
    plt.savefig('gl_min_start.png')
