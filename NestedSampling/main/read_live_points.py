import numpy



filename='results/ganymede/update_nlive/updated_nlive_39.npz'
live_points=numpy.load(filename)

live_u=live_points['arr_0']
live_v=live_points['arr_1']
live_logl=live_points['arr_2']



print(type(live_logl))
