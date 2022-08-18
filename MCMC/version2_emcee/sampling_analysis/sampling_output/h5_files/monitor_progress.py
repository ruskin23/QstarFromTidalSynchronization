from fileinput import filename
import emcee
import numpy

filename="system_10031409.h5"
reader=emcee.backends.HDFBackend(filename,read_only=True)

flatchain=reader.get_chain(flat=True)
last_sample=reader.get_last_sample()

print(flatchain.shape)