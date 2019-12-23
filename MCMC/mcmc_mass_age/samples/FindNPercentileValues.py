import numpy
from NominalValue import Nominal


system='35'
percentile=0.2

parameterKey=['mass','age','feh']
sampleFile = 'MassAgeFehSamples_'+system+'.txt'

NPercentileValues=Nominal(sampleFile,parameterKey,percentile)
print(NPercentileValues())

