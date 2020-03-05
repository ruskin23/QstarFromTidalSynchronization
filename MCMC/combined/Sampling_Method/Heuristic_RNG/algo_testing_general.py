import numpy
import matplotlib.pyplot as plt
from utils import cummulative_distribution
from scipy.stats.stats import pearsonr

class HeuristicRNG:


    def _norm(self,x,mu,sigma):
        return (1/sigma*numpy.sqrt(2*3.14))*numpy.exp(-0.5*((x-mu)/(sigma))**2)

    def _sample(self,mu,sigma):

        return numpy.random.normal(loc=mu,scale=sigma)

    def find_correlations(self):

        R=numpy.zeros(self.desired_covariance.shape)
        l=int(numpy.sqrt(self.desired_covariance.size))

        for i in range(l):
            for j in range(l):
                R[i,j]=self.desired_covariance[i,j]/(numpy.sqrt(self.desired_covariance[i,i]*self.desired_covariance[j,j]))

        return R

    def find_phi_covariance(self):

        C_phi=numpy.zeros([3,3])

        for i in range(3):
            for j in range(3):
                C_phi[i,j]=self.desired_covariance[i,j]

        return C_phi

    def RMSE(self,
             samples=None):

        if samples is None:
            c_qm=numpy.corrcoef(self.phi_matrix[0],self.logQ_samples)[0][1]
            c_qa=numpy.corrcoef(self.phi_matrix[1],self.logQ_samples)[0][1]
            c_qz=numpy.corrcoef(self.phi_matrix[2],self.logQ_samples)[0][1]

        else:
            c_qm=numpy.corrcoef(self.phi_matrix[0],samples)[0][1]
            c_qa=numpy.corrcoef(self.phi_matrix[1],samples)[0][1]
            c_qz=numpy.corrcoef(self.phi_matrix[2],samples)[0][1]

        R=(1/numpy.sqrt(6))*numpy.sqrt((c_qm-self.desired_correlation[3,0])**2 +
                                       (c_qa-self.desired_correlation[3,1])**2 +
                                       (c_qz-self.desired_correlation[3,2])**2)

        return R



    def sample_phi(self):


        N=0
        mean_vector=[]
        for key in self.phi_keys:
            mean_vector=numpy.append(mean_vector,self.parameters[key]['mean'])

        phi_covariance_inverse=numpy.linalg.inv(self.covariance_phi)

        values=[]
        probabilities=[]
        N=0
        with open(self.samples_file,'r') as f:
            next(f)
            for lines in f:
                x=lines.split()
                multiplicity=float(x[3])
                phi_vector=numpy.array([float(x[k]) for k in range(3)])
                arg=0.5*numpy.matmul(numpy.transpose(phi_vector-mean_vector),numpy.matmul(phi_covariance_inverse,phi_vector-mean_vector))
                P=multiplicity*numpy.exp(-arg)

                N=N+P
                values=numpy.append(values,phi_vector)
                probabilities=numpy.append(probabilities,P)


        probabilities=probabilities/N
        values=numpy.reshape(values,(len(values)//3,3))

        U=numpy.random.random(1)[0]

        for k in range(self.constants['sample_size']):
            for i,p in enumerate(probabilities):
                if p>U:
                    self.phi_matrix=numpy.append(self.phi_matrix,values[i])
                    break
                else:U=U-p

        self.phi_matrix=numpy.transpose(numpy.reshape(self.phi_matrix,(self.constants['sample_size'],3)))


    def sample_logQ(self):

        logQ_mean=self.parameters['logQ']['mean']
        logQ_sigma=self.parameters['logQ']['sigma']

        for k in range(self.constants['sample_size']):
            self.logQ_samples=numpy.append(self.logQ_samples,self._sample(logQ_mean,logQ_sigma))



    def Algo(self,RMSE0):

        M=self.constants['sample_size']
        T=self.constants['T']
        dT=self.constants['dT']
        alpha=self.constants['alpha']
        beta=self.constants['beta']
        gamma=self.constants['gamma']

        RMSE_i=RMSE0
        logQ_best=numpy.copy(self.logQ_samples)
        logQ_best_prev=numpy.copy(self.logQ_samples)

        N=0
        i=0

        while i<alpha:
            n_total=0
            n_accepted=0
            for k in range(beta*M):
                N=N+1
                n_total=n_total+1

                i1=numpy.random.randint(0,M)
                i2=numpy.random.randint(0,M)

                self.logQ_samples[i1],self.logQ_samples[i2]=self.logQ_samples[i2],self.logQ_samples[i1]

                RMSE_new=self.RMSE()

                delta=RMSE_i-RMSE_new

                with open(self.analysis_file,'a') as f:
                    f.write(repr(i)+'\t'+
                            repr(N)+'\t'+
                            repr(RMSE_i)+'\t'+
                            repr(RMSE_new)+'\t'+
                            repr(delta)+'\t')

                if delta>0 or delta==0:
                    RMSE_i=RMSE_new
                    n_accepted=n_accepted+1
                    with open(self.analysis_file,'a') as f:
                        f.write('accepted'+'\n')
                else:
                    U=numpy.random.random(1)[0]
                    P=numpy.exp(delta/T)
                    if P>U:
                        RMSE_i=RMSE_new
                        n_accepted=n_accepted+1
                        with open(self.analysis_file,'a') as f:
                            f.write('C_accepted'+'\n')
                    else:
                        self.logQ_samples[i2],self.logQ_samples[i1]=self.logQ_samples[i1],self.logQ_samples[i2]
                        with open(self.analysis_file,'a') as f:
                            f.write('rejected'+'\n')

            if self.RMSE()<self.RMSE(samples=logQ_best):
                logQ_best=numpy.copy(self.logQ_samples)

            if self.RMSE(samples=logQ_best)<self.RMSE(logQ_best_prev):
                i=0

            r=n_accepted/n_total
            if r<gamma:
                i=i+1

            logQ_best_prev=numpy.copy(logQ_best)

            T=dT*T


    def __init__(self,
                 system,
                 parameters,
                 desired_covariance,
                 constants):

        self.system=system
        self.parameters=parameters
        self.constants=constants
        self.desired_covariance=desired_covariance

        self.desired_correlation=self.find_correlations()

        self.phi_keys=['mass','age','feh']
        self.covariance_phi=self.find_phi_covariance()

        self.samples_file='../../../mcmc_mass_age/samples/updated_samples/MassAgeFehSamples_'+self.system+'.txt'

        self.phi_matrix=[]
        self.logQ_samples=[]

        self.analysis_file='/mnt/md0/ruskin/QstarFromTidalSynchronization/sampling_algo/analysis.txt'



    def __call__(self):

        self.sample_phi()
        self.sample_logQ()

        RMSE0=self.RMSE()

        with open(self.analysis_file,'w') as f:
            f.write('i'+'\t'+
                    'N'+'\t'+
                    'RMSE_i'+'\t'+
                    'RMSE_new'+'\t'+
                    'delta'+'\t'+
                    'a or r'+'\n'
                    )

        self.Algo(RMSE0)

        print(self.RMSE())

        return self.phi_matrix,self.logQ_samples




if __name__=='__main__':


    system='137'
    accepted_filename='../../SAVED_CHAINS/ganymede/MCMC_'+system+'/accepted_parameters_1.txt'

    with open(accepted_filename,'r') as f:
        for lines in f:
            x=lines.split()

    parameters=dict(logQ=dict(mean=float(x[4]),
                              sigma=0.5),
                    mass=dict(mean=float(x[5]),
                              sigma=1.0),
                    age=dict(mean=float(x[6]),
                             sigma=1.0),
                    feh=dict(mean=float(x[7]),
                             sigma=0.3)
                    )

    desired_covariance=numpy.array([[0.00577367,-0.064245,0.00930119,0.01271759],
                                    [-0.064245,4.68224828,-0.12771926,0.08551816],
                                    [0.00930119,-0.12771926,0.0307736,0.00915663],
                                    [0.01271759,0.08551816,0.00915663 ,0.08012553]])

    constants=dict(sample_size=100,
                   T=0.01,
                   dT=0.9,
                   alpha=10,
                   beta=5,
                   gamma=0.01)

    O=HeuristicRNG(system,
                        parameters,
                        desired_covariance,
                        constants)
    PHI,LOGQ=O()
