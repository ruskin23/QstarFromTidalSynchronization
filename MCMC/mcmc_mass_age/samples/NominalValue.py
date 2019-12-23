from CummulativeDistribution import CummulativeDistribution
import argparse
import numpy

class Nominal:

    def NPercentileValue(self):

        Values={}
        Distribution=CummulativeDistribution(self.SampleFile,self.PercentileValue)
        for parameter in self.parameterKey:
            Values[parameter]=Distribution(parameter,option='PercentileValue')
        return Values

    def CorrectNominal(self,
                       Values):

        SquaredValues=[]
        parameterDict={}

        for parameter in self.parameterKey:
            parameterDict[parameter]=[]

        with open(self.SampleFile,'r') as f:
            for lines in f:
                x=lines.split()
                SquaredSum=0
                for parameter in self.parameterKey:
                    parameterValue=float(x[self.parameterKey.index(parameter)+1])
                    SquaredSum=SquaredSum+(((Values[parameter][0]-parameterValue)/Values[parameter][1])**2)
                    parameterDict[parameter].append(parameterValue)
                SquaredValues.append(numpy.sqrt(SquaredSum))

        MinSquaredIndex=SquaredValues.index(min(SquaredValues))

        NominalValues=[]
        for parameter in self.parameterKey:
            NominalValues.append(parameterDict[parameter][MinSquaredIndex])

        return NominalValues

    def __call__(self):

        StandardMean=self.NPercentileValue()
        FinalMean=self.CorrectNominal(StandardMean)

        return FinalMean

    def __init__(self,
                 SampleFile,
                 parameterKey,
                 PercentileValue):

        self.SampleFile=SampleFile
        self.parameterKey=parameterKey
        self.PercentileValue=PercentileValue
