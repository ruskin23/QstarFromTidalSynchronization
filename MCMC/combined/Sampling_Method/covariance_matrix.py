import numpy


k=1148
d=7

X=numpy.zeros([d,k])

parameter_file='../SAVED_CHAINS/ganymede/AcceptedParameters.txt'




age=[]
logQ=[]
feh=[]
with open(parameter_file,'r') as f:
    for lines in f:
        x=lines.split()
        age=numpy.append(age,float(x[6]))
        logQ=numpy.append(logQ,float(x[4]))
        feh=numpy.append(feh,float(x[7]))


print(numpy.cov(numpy.stack((age,logQ),axis=0)))
C=numpy.cov(numpy.stack((age,logQ),axis=0))
R=numpy.zeros([2,2])
for i in range(2):
    for j in range(2):
        R[i,j]=C[i,j]/(numpy.sqrt(C[i,i]*C[j,j]))


print(R)

"""
with open(parameter_file,'r') as f:
    for i,lines in enumerate(f):
        x=lines.split()
        parameters=numpy.array([float(x[k]) for k in range(1,8)])
        X[0:d,i]=parameters

X_mean=numpy.zeros([d,1])

for i in range(0,d-1):
    X_mean[i,0]=numpy.mean(X[i,0:k-1])

X1=X-X_mean
X2=numpy.transpose(X1)

Sigma=(1/(k-1))*numpy.dot(X1,X2)
print(numpy.dot(Sigma,numpy.transpose(Sigma)))

Samples=numpy.random.multivariate_normal(numpy.transpose(X_mean)[0],Sigma)
#print(Samples)
"""
