import numpy
from NominalValue import Nominal


AllSystems=['1', '2', '4', '5', '7', '8', '10', '12', '13', '14', '15', '17', '18', '20', '23', '25', '26', '27', '28', '29', '30', '31', '32', '33', '35', '36', '37', '38', '39', '40', '41', '43', '44',   '47', '48', '49', '50', '51', '52', '54', '55', '56', '57', '59', '61', '62', '63', '65', '67', '68', '69', '70', '72', '73', '75', '76', '77', '78', '79', '80', '81', '82', '83', '84', '85', '86', '88', '90',  '91', '92', '93', '94', '95', '96', '99', '100', '101', '103', '104', '105', '106', '107', '108', '109', '110', '111', '112', '113', '114', '115', '117', '118', '119', '120', '121', '123', '124', '125', '126',  '128', '129', '130', '131', '132', '133', '134', '137', '138', '139', '140', '141', '142']

parameterKey=['mass','age','feh']

CatalogFile='catalog_KIC.txt'
with open(CatalogFile,'r') as f:
    for i,lines in enumerate(f):
        if i==0:
            x=lines.split()
            print(x)
            Header='\t'.join(x)
            break

UpdatedHead=Header+'\t'+'\t'.join(parameterKey)+'\n'
UpdatedCataogFile='catalog_l0.4.txt'

with open(UpdatedCataogFile,'w',1) as f:
    f.write(UpdatedHead)

for system in AllSystems:

    print('\nSystem = ',system)
    SampleFile='MassAgeFehSamples_'+system+'.txt'
    N=Nominal(SampleFile,parameterKey,0.5)
    NominalValuesObtained=N()
    print('Values Obtained = ',NominalValuesObtained)

    N=[]
    for v in NominalValuesObtained:
        N.append(repr(v))
    NominalWriteUp='\t'.join(N)

    with open(CatalogFile,'r') as f:
        next(f)
        for lines in f:
            x=lines.split()
            atSystem=x[0]
            if system==atSystem:
                x=lines.split()
                Values='\t'.join(x)
                break

    UpdatedParameters=Values+'\t'+NominalWriteUp+'\n'
    with open(UpdatedCataogFile,'a',1) as fN:
        fN.write(UpdatedParameters)
