import numpy
import itertools

def _Norm(value,loc=0.0,sigma=1.0):
    arg=(value-loc)/sigma
    return numpy.exp(-(arg**2)/2)


def _multivariate_gaussian(x_vector,x_mean,y_vector,y_mean,sigma_xy):
    arg=numpy.matmul(numpy.transpose(x_vector-x_mean),numpy.matmul(sigma_xy,(y_vector-y_mean)))
    return numpy.exp(-0.5*arg)


def _cummulative_distribution(samples,mulitplicity=None):

        if mulitplicity is None: sorted_samples=sorted([(x, len(list(y))) for x, y in itertools.groupby(samples)], key=lambda tup: tup[0])
        else: sorted_samples=sorted(zip(samples,mulitplicity), key=lambda tup: tup[0])

        sample_values=[]
        probability=[]
        p_sum=0.0

        for s in sorted_samples:
            sample_values=numpy.append(sample_values,s[0])
            p_sum=p_sum+s[1]
            probability=numpy.append(probability,p_sum)

        return list(zip(sample_values,probability/max(probability)))


def _get_filename(system,cluster,instance):
#returns the filename with directory of accepted parameter file
    return '../'+cluster+'/MCMC_'+system+'/accepted_parameters_'+str(instance)+'.txt'


def _write_on_file(line,filename,option):
#write line on given filename
#option:'w' to write new file
#       'a' to append
    with open(filename,option) as f:
        f.write(line)


def _fill_parameters(chain_filename):
#creates a file filled_parameters.txt which repeats the missing parameters
#in accepted parameter file

    filled_filename='filled_parameters.txt'
    with open(chain_filename,'r') as f:
        next(f)
        counter=1
        for i,params in enumerate(f):
            x=params.split()
            if i==0:
                # _write_on_file(params,filled_filename,'w')
                # saved_state=x[1:-1]
                # counter=counter+int(x[0])
                # continue
                if x[0]!='1':raise ValueError
                else:
                    _write_on_file(params,filled_filename,'w')
                    saved_state=x[1:-1]
                    counter=counter+1
                    continue
            else:
                iteration=int(x[0])
                if iteration==counter:
                    _write_on_file(params,filled_filename,'a')
                    saved_state=x[1:-1]
                    counter=counter+1
                    continue
                else:
                    difference=iteration-counter
                    for _ in range(difference):
                        line='\t'.join([str(counter)]+saved_state+['\n'])
                        _write_on_file(line,filled_filename,'a')
                        counter=counter+1
                    _write_on_file(params,filled_filename,'a')
                    saved_state=x[1:-1]
                    counter=iteration+1


    return filled_filename


