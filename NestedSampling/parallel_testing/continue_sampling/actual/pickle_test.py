import dynesty
import dill





with open('initial_sampling_saved.dill','rb') as f:
    dsampler=dill.load(f)

print(dsampler.sampler.live_u)
