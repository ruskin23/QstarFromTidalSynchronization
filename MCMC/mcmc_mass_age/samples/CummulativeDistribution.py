import numpy as np
import matplotlib.pyplot as plt
import sys


class CummulativeDistribution():

    def Historgram(self,
                   values):

        n_bins=1000

        fig, ax = plt.subplots(figsize=(8, 4))
        n, bins, patches = ax.hist(values, bins=n_bins, density=True, cumulative=True, histtype='step')
        return n,bins

    def __call__(self,
                 parameter=None,
                 option=None):

        values=[]
        ParameterKey=['mass','age','feh']
        Index=ParameterKey.index(parameter)+1

        with open(self.SampleFile,'r') as f:
            for lines in f:
                x=lines.split()
                values.append(float(x[Index]))

        values.sort()
        n,bins=self.Historgram(values)

        if option=='PercentileValue':
            D=dict(zip(n,bins))
            for key,value in D.items():
                if key>self.percentile:
                    Mean=value
                    break


        Std=np.std(values)
        return Mean,Std

        if option=='Plot':
            ax.plot(bins, y, 'k--', linewidth=1.5)
            ax.grid(True)
            ax.set_title('Cumulative step histograms')
            ax.set_xlabel('age')
            ax.set_ylabel('Likelihood of occurrence')
            plt.show()


    def __init__(self,
                 SampleFile,
                 percentile=None):

        self.SampleFile=SampleFile
        self.percentile=percentile
"""

if __name__ == '__main__':


    #AllSystems=['1', '2', '4', '5', '7', '8', '10', '12', '13', '14', '15', '17', '18', '20', '23', '25', '26', '27', '28', '29', '30', '31', '32', '33', '35', '36', '37', '38', '39', '40', '41', '43', '44', '47', '48', '49', '50', '51', '52', '54', '55', '56', '57', '59', '61', '62', '63', '65', '67', '68', '69', '70', '72', '73', '75', '76', '77', '78', '79', '80', '81', '82', '83', '84', '85', '86', '88', '90', '91', '92', '93', '94', '95', '96', '99', '100', '101', '103', '104', '105', '106', '107', '108', '109', '110', '111', '112', '113', '114', '115', '117', '118', '119', '120', '121', '123', '124', '125', '126', '128', '129', '130', '131', '132', '133', '134', '137', '138', '139', '140', '141', '142']

    AllSystems=['2']
    for system in AllSystems:
        print('\nSystem = ',system)
        SampleFile='MassAgeFehSamples_'+system+'.txt'
        parameter='age'
        distribution=CummulativeDistribution(SampleFile,0.5)
        D={parameter:distribution(parameter=parameter,option='PercentileValue')}

    print(D['age'][0])


"""
