import numpy
from utils import _get_filename, _fill_parameters, _get_chain, adjust_chain
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec



def z_value(chain_a,chain_b):

    mean_a=numpy.mean(chain_a)
    var_a=numpy.var(chain_a)

    mean_b=numpy.mean(chain_b)
    var_b=numpy.var(chain_b)

    return (mean_a-mean_b)/(numpy.sqrt(var_a+var_b))

def break_chain(chain):

    #break chain into initial 20%, in between values and last 50%
    initial_part,discard_part,last_part=numpy.split(chain,[int(len(chain)*0.1),int(len(chain)*0.5)])

    div=[int(len(last_part)*0.05*k) for k in range(1,20)]
    #divide initial part into 10 chains
    last_arrays=numpy.split(last_part,div)

    return div,initial_part,last_arrays

if __name__=='__main__':

    #systems=['8','43','36','109','70','47','86','88','93','123','95','106','79','84','25','12','50','28','13']
    systems=['1', '8', '12', '13', '17', '20', '25', '28', '32', '36', '39', '43', '44', '47', '48', '50', '54', '56', '67', '70', '73', '76', '79', '80', '81', '83', '84', '85', '86', '88', '92', '93', '94', '95', '96', '106', '109', '120', '123', '126', '137']


    plt.figure(figsize=(15,10))
    gs1=gridspec.GridSpec(6,7)
    gs1.update(wspace=0.0, hspace=0.0)

    for s in range(41):

        print('System = ',systems[s])
        INITIAL=numpy.zeros(1)
        FINAL_SEGMENTS=[numpy.zeros(1) for k in range(20)]
        for c in ['ganymede','stampede']:
            for i in range(1,6):
                chain_filename=_get_filename(systems[s],c,i)
                filled_name=_fill_parameters(chain_filename)
                chain=_get_chain('logQ',filled_name)
                chain=adjust_chain(systems[s],chain)
                div,initial_part,last_arrays=break_chain(chain)
                INITIAL=numpy.concatenate((INITIAL,initial_part),axis=None)
                if INITIAL[0]==0:INITIAL=INITIAL[1:]
                for i,f in enumerate(last_arrays):
                    FINAL_SEGMENTS[i]=numpy.concatenate((FINAL_SEGMENTS[i],f),axis=None)
                    if FINAL_SEGMENTS[i][0]==0:FINAL_SEGMENTS[i]=FINAL_SEGMENTS[i][1:0]
        
        Z=[]
        IT=[]
        for i,f in enumerate(FINAL_SEGMENTS):
            Z.append(z_value(INITIAL,f))
            IT.append(i)

        plt.axis('on')
        ax1=plt.subplot(gs1[s])
        ax1.scatter(IT,Z)
        ax1.hlines(2.0,min(IT),max(IT),linestyles='dashed')
        ax1.hlines(-2.0,min(IT),max(IT),linestyles='dashed')
        ax1.tick_params(axis='x', labelsize= 5)
        ax1.tick_params(axis='y', labelsize= 5)
        ax1.label_outer()

    plt.savefig('all_pdf.png')

    # for system in systems:
    #     print('System = ',system)
    #     INITIAL=numpy.zeros(1)
    #     FINAL_SEGMENTS=[numpy.zeros(1) for k in range(20)]
    #     for c in ['ganymede','stampede']:
    #         for i in range(1,6):
    #             chain_filename=_get_filename(system,c,i)
    #             filled_name=_fill_parameters(chain_filename)
    #             chain=_get_chain('logQ',filled_name)
    #             chain=adjust_chain(systems,chain)
    #             div,initial_part,last_arrays=break_chain(chain)
    #             INITIAL=numpy.concatenate((INITIAL,initial_part),axis=None)
    #             if INITIAL[0]==0:INITIAL=INITIAL[1:]
    #             for i,f in enumerate(last_arrays):
    #                 FINAL_SEGMENTS[i]=numpy.concatenate((FINAL_SEGMENTS[i],f),axis=None)
    #                 if FINAL_SEGMENTS[i][0]==0:FINAL_SEGMENTS[i]=FINAL_SEGMENTS[i][1:0]
        
    #     Z=[]
    #     IT=[]
    #     for i,f in enumerate(FINAL_SEGMENTS):
    #         Z.append(z_value(INITIAL,f))
    #         IT.append(i)
    #     plt.scatter(IT,Z)
    #     plt.hlines(2.0,min(IT),max(IT),linestyles='dashed')
    #     plt.hlines(-2.0,min(IT),max(IT),linestyles='dashed')
    #     plt.savefig(f'plots/z_test_new/System_{system}.png')
    #     plt.close()
