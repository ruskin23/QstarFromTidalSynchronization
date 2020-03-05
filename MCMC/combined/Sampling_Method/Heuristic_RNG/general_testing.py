import numpy
import matplotlib.pyplot as plt
from utils import cummulative_distribution
from covariance_matrix import Covariance
from scipy.stats.stats import pearsonr
import time
import sys

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
             theta,
             samples=None):

        r=self.current_theta_index+1
        c=3+self.current_theta_index

        R=numpy.zeros([r,c])

        for i in range(r):
            for j in range(3+i):
                R[i,j]=numpy.corrcoef(self.Y_matrix[j],self.Y_matrix[3+i])-self.desired_correlation[j,3+i                           ]

        RT=numpy.transpose(R)

        rmse=0
        for i in range(r):
            rmse=Rmse+numpy.dot(R[i,:],RT[:,i])

        return rmse


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


        for k in range(self.constants['sample_size']):
            U=numpy.random.random(1)[0]
            for i,p in enumerate(probabilities):
                if p>U:
                    self.phi_matrix=numpy.append(self.phi_matrix,values[i])
                    break
                else:U=U-p

        self.phi_matrix=numpy.transpose(numpy.reshape(self.phi_matrix,(self.constants['sample_size'],3)))


        print('Mass_age =',numpy.corrcoef(self.phi_matrix[0],self.phi_matrix[1])[0][1]-self.desired_correlation[2,1])

    def sample_theta(self):

        theta_mean=self.parameters[self.current_theta]['mean']
        theta_sigma=self.parameters[self.current_theta]['sigma']

        for k in range(self.constants['sample_size']):
            self.theta_samples=numpy.append(self.theta_samples,self._sample(theta_mean,theta_sigma))



    def Algo(self,RMSE0):

        M=self.constants['sample_size']
        T=self.constants['T']
        dT=self.constants['dT']
        alpha=self.constants['alpha']
        beta=self.constants['beta']
        gamma=self.constants['gamma']

        RMSE_i=RMSE0

        yk=3+self.current_theta_index

        theta_best=numpy.copy(self.Y_matrix[yk])
        theta_best_prev=numpy.copy(theta_best)

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

                self.Y_matrix[yk,i1],self.Y_matrix[yk,i2]=self.Y_matrix[yk,i2],self.Y_matrix[yk,i1]

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
                        self.Y_matrix[yk,i2],self.Y_matrix[yk,i1]=self.Y_matrix[yk,i1],self.Y_matrix[yk,i2]
                        with open(self.analysis_file,'a') as f:
                            f.write('rejected'+'\n')


            if self.RMSE()<self.RMSE(samples=theta_best):
                theta_best=numpy.copy(self.logQ_samples)

            if self.RMSE(samples=theta_best)<self.RMSE(samples=logQ_best_prev):
                i=0

            r=n_accepted/n_total
            if r<gamma:
                i=i+1

            theta_best_prev=numpy.copy(logQ_best)

            T=dT*T


    def __init__(self,
                 instance,
                 system,
                 parameters,
                 constants):


        self.instance=instance
        self.system=system
        self.parameters=parameters
        self.constants=constants

        C=Covariance(system)

        self.desired_covariance=C.Calculate('Covariance')
        self.desired_correlation=C.Calculate('Correlation')

        self.phi_keys=['mass','age','feh']
        self.covariance_phi=self.find_phi_covariance()

        self.theta_keys=['Porb','eccentricity','Wdisk','logQ']
        self.current_theta=None
        self.current_theta_index=None

        self.samples_file='../../../mcmc_mass_age/samples/updated_samples/MassAgeFehSamples_'+self.system+'.txt'

        self.phi_matrix=[]
        self.theta_samples=[]

        self.analysis_file='/mnt/md0/ruskin/QstarFromTidalSynchronization/sampling_algo/analysis_'+self.instance+'.txt'

        with open(self.analysis_file,'w') as f:
            f.write('i'+'\t'+
                    'N'+'\t'+
                    'RMSE_i'+'\t'+
                    'RMSE_new'+'\t'+
                    'delta'+'\t'+
                    'a or r'+'\n'
                    )


    def __call__(self):

        start_time = time.time()

        self.sample_phi()

        self.Y_matrix=numpy.copy(self.phi_matrix)

        for keys in self.theta_keys:

            print('Calculating for theta = ',keys)
            self.theta_samples=[]

            self.current_theta=keys
            self.current_theta_index=self.theta.index(keys)

            self.sample_theta()

            self.Y_matrix=numpy.vstack(self.Y_matrix,self.theta_samples)

            RMSE0=self.RMSE()

            self.Algo(RMSE0)

            print('Phi = ',self.phi_matrix)
            print('R = ',self.desired_correlation)
            print(self.RMSE())
            print("--- %s seconds ---" % (time.time() - start_time))

            C=numpy.cov(self.Y_matrix)
            print(C)


            R=numpy.zeros([7,7])
            for i in len(7):
                for j in len(7):
                    R[i,j]=C[i,j]/numpy.sqrt(C[i,i]*C[j,j])

            print(R)


        return self.Y_matrix,R




if __name__=='__main__':

    instance=sys.argv[1]
    system='137'
    accepted_filename='../../SAVED_CHAINS/ganymede/MCMC_'+system+'/accepted_parameters_1.txt'

    with open(accepted_filename,'r') as f:
        for lines in f:
            x=lines.split()

    parameters=dict(theta=dict(mean=float(x[4]),
                              sigma=0.5),
                    mass=dict(mean=float(x[5]),
                              sigma=1.0),
                    age=dict(mean=float(x[6]),
                             sigma=1.0),
                    feh=dict(mean=float(x[7]),
                             sigma=0.3)
                    )

    constants=dict(sample_size=1000,
                   T=0.01,
                   dT=0.9,
                   alpha=5,
                   beta=5,
                   gamma=0.01)

    O=HeuristicRNG(instance,
                    system,
                    parameters,
                    constants)
    PHI,LOGQ=O()
