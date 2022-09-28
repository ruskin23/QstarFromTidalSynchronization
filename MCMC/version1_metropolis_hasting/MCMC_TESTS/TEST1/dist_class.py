class CummulativeDistribution:

    def __init__(self,
                 filename,
                 param):
        self.filename=filename
        self.index=param
        self.data=[]
        self.normalization=-1

    def generate_data(self):
        multiplicity=0
        with open(sample_file_name,'r') as f:
            next(f)
            for i,lines in enumerate(f):
                x=lines.split()
                current_state=[float(x[k]) for k in range(1,5)]
                if i==0:previous_state=current_state
                if previous_state==current_state:
                    multiplicity=multiplicity+1
                else:
                    previous_state.append(multiplicity)
                    self.data.append(previous_state)
                    previous_state=current_state
                    multiplicity=1
                self.normalization=self.normalization+1

    def sort_data(self):
        sorted_DATA=sorted(self.data,key = lambda d: d[self.index])
        return sorted_DATA

    def __call__(self):
        self.generate_data()

        value=[]
        m=[]
        n=0

        for d in self.sort_data():
            value.append(d[param])
            n=n+d[-1]
            m.append(n/self.normalization)

        return value,m






