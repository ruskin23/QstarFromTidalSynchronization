import numpy


system='39'

filename='results/ganymede/nlive_39.npz'
live_points=numpy.load(filename)

live_u=live_points['arr_0']
live_v=live_points['arr_1']
live_logl=live_points['arr_2']

for i,k in enumerate(live_logl):
    if k==0:
        print(live_logl[i])
        print(i)
        break
print(live_logl[80])