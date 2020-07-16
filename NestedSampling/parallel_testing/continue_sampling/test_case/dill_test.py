import pickle
import dill

with open('sample_initial.dill','rb') as f:
    S=dill.load(f)

print(len(S.results['samples']))
print(len(S.results['samples'][0]))
