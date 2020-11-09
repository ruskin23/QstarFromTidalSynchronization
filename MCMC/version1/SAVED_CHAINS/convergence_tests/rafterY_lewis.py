import pymc
from pymc import raftery_lewis
from utils import _get_filename, _fill_parameters, _get_chain
import numpy

systems=['8','43','36','109','70','47','86','88','93','123','95','106','79','84','25','12','50','28','13']

for system in systems:
    print('\n\n System=',system)
    CHAIN=numpy.zeros(1)
    for c in ['ganymede','stampede']:
        for i in range(1,6):
            chain_filename=_get_filename(system,c,i)
            filled_name=_fill_parameters(chain_filename)
            chain=_get_chain('logQ',filled_name)
            CHAIN=numpy.concatenate((CHAIN,chain),axis=None)

    raftery_lewis(CHAIN,q=0.95, r=0.01)
