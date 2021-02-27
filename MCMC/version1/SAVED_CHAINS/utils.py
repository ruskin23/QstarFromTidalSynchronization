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


def _get_chain(parameter,chain_filename):

    parameters=['Porb','eccentricity','Wdisk','logQ','primary_mass','age','feh']
    if parameter=='Spin':param_idx=15
    elif parameter in parameters:param_idx=parameters.index(parameter)+1
    else:
        print('Parameters names can only be: ',parameters+['Spin'])
        raise ValueError

    CHAIN=[]
    with open(chain_filename,'r') as f:
        next(f)
        for lines in f:
            x=lines.split()
            try:
                CHAIN=numpy.append(CHAIN,float(x[param_idx]))
            except:
                print(chain_filename)
                print(x)
                raise Exception
    return CHAIN


def _get_filename(system,cluster,instance):
#returns the filename with directory of accepted parameter file
    f = '../'+cluster+'/MCMC_'+system+'/accepted_parameters_'+str(instance)+'.txt'
    return f

def _write_on_file(line,filename,option):
#write line on given filename
#option:'w' to write new file
#       'a' to append
    with open(filename,option) as f:
        f.write(line)


def _fill_parameters(chain_filename,skip_first=False):
#creates a file filled_parameters.txt which repeats the missing parameters
#in accepted parameter file

    filled_filename='filled_parameters.txt'
    
    with open(chain_filename,'r') as f:
        next(f)
        if skip_first==True:next(f)
        counter=1
        for i,params in enumerate(f):
            x=params.split()
            if i==0:
                _write_on_file(params,filled_filename,'w')
                saved_state=x[1:-1]
                counter=counter+1
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


def _get_filled_chain(system,parameter='logQ'):

    CHAIN=numpy.zeros(1)
    
    saved_directory='/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/version1/SAVED_CHAINS'
    for c in ['ganymede','stampede']:
        for i in range(1,6):
            chain_filename=f'{saved_directory}/{c}/MCMC_{system}/accepted_parameters_{i}.txt'
            if chain_filename!=f'{saved_directory}/ganymede/MCMC_{system}/accepted_parameters_1.txt':skip_first=True
            else:skip_first=False
            filled_file=_fill_parameters(chain_filename,skip_first=skip_first)
            chain=_get_chain(parameter,filled_file)
            CHAIN=numpy.concatenate((CHAIN,chain),axis=None)
    
    return CHAIN[1:]

def adjust_chain(system,chain):

    delete_q=[]
    saved_directory='/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/version1/SAVED_CHAINS'

    with open(f'{saved_directory}/deletion_checks/final_results.txt','r') as f:
        for lines in f:
            x=lines.split()
            s=x[0].split('_')[0]
            if s==system:
                if x[6]=='True':
                    delete_q.append(float(x[2]))
                else:delete_q.append(8.353237442958045)

    delete_idx=numpy.zeros(1)
    for q in delete_q:
        idx=numpy.where(chain==q)
        delete_idx=numpy.concatenate((delete_idx,idx[0]),axis=None)
    delete_idx=delete_idx.astype(int)
    chain=numpy.delete(chain,delete_idx)
    return chain


def _adjust_mulitplcity(system,samples=None,dump=False):

    # if samples is None:samples=_get_filled_chain(system)
    # sorted_samples=sorted([(x, len(list(y))) for x, y in itertools.groupby(samples)], key=lambda tup: -tup[1])

    # with open(f'deletion_checks/system_{system}.txt','w') as f:
    #     for ss in sorted_samples:
    #         f.write(repr(ss[0])+'\t'+repr(ss[1])+'\n')

    saved_directory='/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/version1/SAVED_CHAINS'    
    for c in ['ganymede','stampede']:
        for i in range(1,6):
            chain_filename=f'{saved_directory}/{c}/MCMC_{system}/accepted_parameters_{i}.txt'
            print(f'{c}_accepted_parameters_{i}.txt')
            filled_file=_fill_parameters(chain_filename)
            chain=_get_chain('logQ',filled_file)
            sorted_samples=sorted([(x, len(list(y))) for x, y in itertools.groupby(chain)], key=lambda tup: -tup[1])
            with open(f'deletion_checks/{system}_{c}_{i}.txt','w') as f:
                for ss in sorted_samples:
                    f.write(repr(ss[0])+'\t'+repr(ss[1])+'\n')
        

