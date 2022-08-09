
from collections import namedtuple
import pickle
from sys import api_version
import numpy
import emcee
# config_dict = dict(num_processes=16)
# config = namedtuple("CONFIG",config_dict.keys())(*config_dict.values())

# print(config.num_processes)


from configargparse import ArgumentParser, DefaultsFormatter


def cmd_parser():
    p = ArgumentParser(description='test',
                       default_config_files=['config.txt'])
    p.add_argument('--a',type=float,default=111,help='test argument')

    return p.parse_args()

def generate_lines_that_equal(string, fp):
    for line in fp:
        x=line.split()
        if x[0] == string:
            yield x[0]

if __name__=='__main__':

    # config=cmd_parser()
    # print(config.a)

    # d=dict(a=1,b=2,c=3)

    
    # a=[]
    # b=(12,)
    # p = tuple(value for key,value in d.items())
    # print(p)
    # for key,values in d.items():
    #     a.append(values)
    # a=tuple(a)
    # print(tuple(
    #     -numpy.inf if i == 0 else numpy.nan
    #     for i in range(10)
    # )) 

    # blobs_dtype = [(name, float) for name, _ in d.items()]
    # print(blobs_dtype)
    # blobs_dtype = numpy.dtype(blobs_dtype)
    # print(blobs_dtype)

    
    # with open('/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/version2_emcee/catalog/filtering/Lurie_binaries_with_p1.txt','r') as f:
    #     for k in generate_lines_that_equal('10264202',f):
    #         KIC=k
    #         break
    
    # reader = emcee.backends.HDFBackend("test_run.h5", read_only=True)
    # flatchain = reader.get_blobs(flat=True)
    # print(flatchain)


    # chain=dict()

    # a=[1,2,3,4]
    # b=[4,5,6,7]
    
    # saved=dict()
    # saved['a']=a
    # saved['b']=b

    # crate=dict()
    # crate['system']=saved


    # with open('test.pickle','wb') as f:
    #     pickle.dump(crate,f,protocol=pickle.HIGHEST_PROTOCOL)

    # with open('test.pickle','rb') as f:
    #     pick=pickle.load(f)

    # print(pick)

    with open('bandwidth.pickle','rb') as f:
        p=pickle.load(f)
    print(p['9971475']['ISJ'])