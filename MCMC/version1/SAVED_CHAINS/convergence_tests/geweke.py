import numpy
from utils import _get_filename, _fill_parameters, _get_chain
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


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
    systems=['85', '73', '76', '96', '92', '81', '80', '36', '93', '83', '84', '94', '32', '79', '106', '123', '50', '47', '39', '56', '126', '54', '109', '44', '48', '17', '70', '8', '12', '88', '67', '20', '95', '25', '57', '137', '120', '86', '43', '28', '13']

    for system in systems:
        print('System = ',system)
        INITIAL=numpy.zeros(1)
        FINAL_SEGMENTS=[numpy.zeros(1) for k in range(20)]
        for c in ['ganymede','stampede']:
            for i in range(1,6):
                chain_filename=_get_filename(system,c,i)
                filled_name=_fill_parameters(chain_filename)
                chain=_get_chain('logQ',filled_name)
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
        plt.scatter(IT,Z)
        plt.hlines(2.0,min(IT),max(IT),linestyles='dashed')
        plt.hlines(-2.0,min(IT),max(IT),linestyles='dashed')
        plt.savefig(f'plots/z_test/System_{system}.png')
        plt.close()
