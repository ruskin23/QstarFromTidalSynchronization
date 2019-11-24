import numpy
import argparse



class TrueNominal():


    def _MeanValues(self):

        with open(self.nominal_file,'r') as f:
            next(f)
            print('defining dictionary for mean value')
            for line in f:
                x=line.split()
                if x[0]==self.system:
                    print('System = ',i)
                    for k,param in enumerate(self.parameter_set):
                        print('Adding for ', param)
                        self.nominal_values[param+'_mean']=float(x[2*k+1])
                        self.nominal_values[param+'_sigma']=float(x[2*k+2])
                    print(self.nominal_values)
                    break

    def CalculateDev(self,
                     parameters):

        dev=0
        for i,param in enumerate(self.parameter_set):
            dev=dev+((parameters[i] - self.nominal_values[param+'_mean'])
                     /
                     self.nominal_values[param+'_sigma']
                     )**2

        return numpy.sqrt(dev)

    def CalculateMin(self):

        self._MeanValues()

        check_dev=[]
        print('Calculating min')
        with open(self.sample_file,'r') as f:
            next(f)
            for k,lines in enumerate(f):
                p=[]
                x=lines.split()
                for i in range(1,len(x)):
                    p.append(float(x[i]))
                print('Parameter Array = ')
                print(p)

                dev=self.CalculateDev(p)

                check_dev.append(dev)

                if k==0:
                    min_dev=dev
                    min_parameters=p
                else:
                    if dev<min_dev:
                        min_dev=dev
                        min_parameters=p
                    else:continue
        print(min_dev)
        print(min(check_dev))
        return min_parameters

    def __init__(self,
                 parameter_set,
                 sample_file,
                 nominal_file,
                 system):

        self.parameter_set=parameter_set

        self.sample_file=sample_file
        self.nominal_file=nominal_file

        self.system=system

        self.nominal_values=dict()


if __name__=='__main__':



    parameter_set=['mass','age','feh']
    nominal_file='nominal_values.txt'

    index=[]
    with open('catalog_KIC.txt','r') as f:
        next(f)
        for lines in f:
            x=lines.split()
            index.append(x[0])




    with open('Corrected_Nominal.txt','w',1) as f:

        head='\t'.join(parameter_set)
        f.write('0'+'\t'+head+'\n')

        for i in index:

            if i in ['111','139']:continue

            sample_file='MassAgeFehSamples_'+i+'.txt'

            correct_mean=TrueNominal(parameter_set,
                                     sample_file,
                                     nominal_file,
                                     i).CalculateMin()


            values=[repr(c) for c in correct_mean]
            val='\t'.join(values)
            f.write(i+'\t'+val+'\n')


