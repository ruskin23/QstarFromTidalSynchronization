import numpy
from NominalValue import Nominal

parameterKey=['mass','age','feh']
sampleFile = 'MassAgeFehSamples_62.txt'
percentile=0.1

NPercentileValues=Nominal(sampleFile,parameterKey,percentile)
print(NPercentileValues()[1])


"""

s=['62', '69', '59', '129', '5', '111', '99', '18', '27', '33', '49', '132', '128', '37', '125', '68']
p=[0.1,0.2,0.3,0.4]

with open('PercentileAges.txt','w',1) as f:
    f.write('system'+'\t'+
            '0.1'+'\t'+
            '0.2'+'\t'+
            '0.3'+'\t'+
            '0.4'+'\n')
    for system in s:
        f.write(system+'\t')
        for percentile in p:
            parameterKey=['mass','age','feh']
            sampleFile = 'MassAgeFehSamples_'+system+'.txt'

            NPercentileValues=Nominal(sampleFile,parameterKey,percentile)
            age=NPercentileValues()[1]
            f.write(repr(age)+'\t')
        f.write('\n')

"""
