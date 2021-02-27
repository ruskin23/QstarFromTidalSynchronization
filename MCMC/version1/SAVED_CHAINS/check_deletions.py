import numpy
import utils
import matplotlib.pyplot as plt

s=['1', '8', '12', '13', '17', '20', '25', '28', '32', '36', '39', '43', '44', '47', '48', '50', '54', '56', '57', '67', '70', '73', '76', '79', '80', '81', '83', '84', '85', '86', '88', '92', '93', '94', '95', '96', '106', '109', '120', '123', '126', '137']

def chain_plot():

    redact_value=8.037172897517499
    it=[]
    q=[]
    with open('ganymede/MCMC_32/accepted_parameters_3.txt','r') as f:
        next(f)
        for lines in f:
            x=lines.split()
            if float(x[4])!=redact_value:
                it.append(int(x[0]))
                q.append(float(x[4]))
    plt.plot(it,q)
    plt.show()

for system in s:
    print(f'\n{system}')
    utils._adjust_mulitplcity(system)

# chain_plot()