
from collections import namedtuple
from sys import api_version
import numpy
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

    d=dict(a=1,b=2,c=3)
    a=[]
    b=(12,)

    p = tuple(value for key,value in d.items())
    print(p)
    for key,values in d.items():
        a.append(values)
    a=tuple(a)

    print(b+a)

    print(tuple(
        -numpy.inf if i == 0 else numpy.nan
        for i in range(10)
    )) 



    with open('/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/version2_emcee/catalog/filtering/Lurie_binaries_with_p1.txt','r') as f:
        for k in generate_lines_that_equal('10264202',f):
            KIC=k
            break

    print(KIC)