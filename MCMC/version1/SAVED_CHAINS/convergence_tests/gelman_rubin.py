#system 57 has one chain with 3 accepted. so need to do something bout it

from utils import _get_filename,_get_chain,_fill_parameters
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

    min_len=min(l_chain)
    min_chain_idx=l_chain.index(min_len)

    CHAINS_R=[]
    for _,c in enumerate(CHAINS):
        r=int(len(c)/len(CHAINS[min_chain_idx]))
        c=c[::r]
        diff=len(c)-min_len
        if diff==0:
            c=c[::20]
            CHAINS_R.append(c)
            continue
        r=int(len(c)/diff)
        ridx=[]
        for i in range(len(c)):
            if i%(r-1)==0:ridx.append(i)
            if len(ridx)==diff:break
        c=numpy.delete(c,ridx)
        c=c[::20]
        CHAINS_R.append(c)
    
    return CHAINS_R


def get_chains(system):
    #returns 10 chains of equal length
    print(f'Calulating for system = {system}')
    CHAINS=[]
    l_chain=[]

    for c in ['ganymede','stampede']:
        for i in range(1,6):
            filename=_get_filename(system,c,i)
            filled_file=_fill_parameters(filename)
            chain=_get_chain('logQ',filled_file)
            # chain=chain[::10]
            CHAINS.append(chain)
            l_chain.append(len(chain))

    C=reduce_chain(CHAINS,l_chain)
    
    # for i in range(10):
    #     CHAINS[i]=CHAINS[i][-min_len:-1]

    return C

if __name__=='__main__':

    with open('gr_statistics.txt','w') as f:
         f.write('System'+'\t'+'N'+'\t'+'B'+'\t'+'W'+'\t'+'R'+'\n')
    #potential_systems=['8','43','36','109','70','47','86','88','93','123','95','106','79','84','25','12','50','28','13']
    potential_systems=['70', '36' ,'43' ,'88' ,'86', '8' ,'109', '93', '95' ,'47', '1' ,'123', '106' ,'13', '79', '48' ,'84', '120', '25' ,'12' ,'96' ,'50' ,'28', '94' ,'67']
    non_systems=['76','39','54','80','81','92','126']
    systems=['85', '73', '76', '96', '92', '81', '80', '36', '93', '83', '84', '94', '32', '79', '106', '123', '50', '47', '39', '56', '126', '54', '109', '44', '48', '17', '70', '8', '12', '88', '67', '20', '95', '25', '137', '120', '86', '43', '28', '13']
    N=[]
    R=[]
    colors=[]
    for i,system in enumerate(systems):
        #print(f'\nSystem = {system}')
        chains=get_chains(system)
        n,b,w,r=gr_test(chains)
        N.append(numpy.log10(n))
        R.append(numpy.log10(r))
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
            N.append(numpy.log10(float(x[1])))
            R.append(numpy.log10(float(x[4])))
            if system in non_systems:colors.append('green')
            else:colors.append('red')
    plt.scatter(N,R,color=colors)
    plt.hlines(numpy.log10(1),min(N),max(N),linestyles='dashed',color='k',label=1.0)
    plt.hlines(numpy.log10(1.1),min(N),max(N),linestyles='dashed',color='b',label=1.1)
    plt.hlines(numpy.log10(1.2),min(N),max(N),linestyles='dashed',color='m',label=1.2)
    plt.legend()
    plt.savefig('thin_20.png')
