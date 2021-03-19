import numpy
import pickle
import matplotlib.pyplot as plt
import utils
from scipy.stats import norm
from scipy import interpolate
from scipy import integrate


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

def pdf_calculator(chain):

    x=numpy.linspace(5,12,10000)
    N=len(chain)
    h=3.5*numpy.power(N,-1.0/3)

    f_h=0
    for i in range(N):
        y=(x-chain[i])/h
        f_h=f_h+norm.pdf(y)

    return f_h/(N*h)


# def main_pdf(system):

#     with open('../complete_chains.pickle','rb') as f:
#         D=pickle.load(f)
#     for system_name,params in D.items():
#         if system_name==system:
#             for name,values in params.items():
#                 if name=='logQ':chain=values

#     pdf=pdf_calculator(chain)

#     return pdf

def single_pdf(system):

    main_chain=numpy.zeros(1)
    pdf=[]
    combined_pdf=numpy.ones(10000)

    for c in ['ganymede','stampede']:

        for i in range(1,6):

            chain_filename=utils._get_filename(system,c,i)
            filled_filename=utils._fill_parameters(chain_filename)
            chain=utils._get_chain('logQ',filled_filename)
            print(f'\ngot chain {chain_filename}')
            chain_fixed=adjust_chain(system,chain)
            print('fixed chain')

            p=pdf_calculator(chain_fixed)
            print('single pdf calculated')
            combined_pdf=combined_pdf*p
            pdf.append(p)
            print('combined pdf calculated')


            if chain_filename!=f'ganymede/MCMC_{system}/accepted_parameters_1.txt':
                start_point=chain_fixed[0]
                count_repeat=numpy.count_nonzero(chain_fixed==start_point)
                chain_fixed=chain_fixed[count_repeat:]
            main_chain=numpy.concatenate((main_chain,chain_fixed),axis=None)
    
    main_chain=main_chain[1:]
    
    main_pdf=pdf_calculator(main_chain)

    cdf=utils._cummulative_distribution(main_chain)
    Z=list(zip(*cdf))
    logQ_values=numpy.array(list(Z[0]))
    logQ_cdf=numpy.array(list(Z[1]))
    main_cdf_interp=interpolate.interp1d(logQ_values,logQ_cdf)

    # return pdf,combined_pdf,main_cdf_interp,logQ_values

    return pdf,combined_pdf,main_pdf,main_cdf_interp,logQ_values

def get_quantile(x,main_cdf_interp,logQ_values):
    
    if x<=min(logQ_values): return 0
    elif x>=max(logQ_values): return 1
    else: return main_cdf_interp(x)

def E_p(main_cdf_interp,c_pdf_interp,c_pdf_N,logQ_values):

    I=lambda x: (1/c_pdf_N)*(2*min(get_quantile(x,main_cdf_interp,logQ_values),1-get_quantile(x,main_cdf_interp,logQ_values)))*c_pdf_interp(x)
    p=integrate.quad(I,5,12)[0]
    return p

if __name__=='__main__':
    s=['1', '8', '12', '13', '17', '20', '25', '28', '32', '36', '39', '43', '44', '47', '48', '50', '54', '56', '67', '70', '73', '76', '79', '80', '81', '83', '84', '85', '86', '88', '92', '93', '94', '95', '96', '106', '109', '120', '123', '126', '137']
    with open('E_p_all.txt','w',1) as f:
        for system in s:
            logQ_array=numpy.linspace(5,12,10000)

            # p_list,c_pdf,main_cdf_interp,logQ_values=single_pdf(system)

            p_list,c_pdf,m_pdf,main_cdf_interp,logQ_values=single_pdf(system)

            c_pdf_interp=interpolate.InterpolatedUnivariateSpline(logQ_array,c_pdf)
            c_pdf_N=c_pdf_interp.integral(5,12)

            main_pdf_interp=interpolate.InterpolatedUnivariateSpline(logQ_array,m_pdf)
            main_pdf_N=main_pdf_interp.integral(5,12)

            p = E_p(main_cdf_interp,c_pdf_interp,c_pdf_N,logQ_values)

            f.write(system+'\t'+repr(p)+'\n')
            plt.plot(logQ_array,m_pdf/main_pdf_N,color='black')
            plt.plot(logQ_array,c_pdf/c_pdf_N,color='red')
            # for p in p_list:;
            #     plt.plot(logQ_array,p/max(p),linestyle='--')
            plt.savefig(f'pdf/system_{system}.png')
            plt.close()

