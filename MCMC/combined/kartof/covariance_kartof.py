import numpy

class Covariance:

    def __init__(self,
                 system):

        self.system=system
        self.parameters_keys=['Porb','e','Wdisk','logQ','mass','age','feh']

        self.parameters=dict()
        for key in self.parameters_keys:
            self.parameters[key]=[]

    def get_samples(self):

        parameter_file='AcceptedParameters.txt'

        with open(parameter_file,'r') as f:
            for lines in f:
                x=lines.split()

                for i, key in enumerate(self.parameters_keys):
                    self.parameters[key]=numpy.append(self.parameters[key],float(x[i+1]))


    def Calculate(self,
                  required):

        self.get_samples()
        parameters_matrix=numpy.array([self.parameters['Porb'],
                                       self.parameters['e'],
                                       self.parameters['Wdisk'],
                                       self.parameters['logQ'],
                                       self.parameters['mass'],
                                       self.parameters['age'],
                                       self.parameters['feh']
                                       ])

        C=numpy.cov(parameters_matrix)
        L=7
        R=numpy.zeros([L,L])
        for i in range(L):
            for j in range(L):
                R[i,j]=C[i,j]/(numpy.sqrt(C[i,i]*C[j,j]))

        if required=='Covariance':return C
        if required=='Correlation':return R
