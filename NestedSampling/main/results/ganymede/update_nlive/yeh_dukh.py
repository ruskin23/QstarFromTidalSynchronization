import numpy

systems = ['39', '54', '76', '92', '126']

for s in systems:
    logl=numpy.zeros(500)

    with open(f'output_{s}.txt','r') as f:
        for lines in f:
            x=lines.split()
            if x[0]=='Running':
                i_min,i_max=int(x[3]),int(x[5])-1
                i=i_min
            if x[0]=='loglike:':
                logl[i]=float(x[1])
                i=i+1


    filename='../nlive_39.npz'
    live_points=numpy.load(filename)
    live_u=live_points['arr_0']
    live_v=live_points['arr_1']
    live_logl=live_points['arr_2']

    for i,l in enumerate(live_logl):
        if l==0:live_logl[i]=logl[i]

    outfile='updated_nlive_39.npz'
    numpy.savez(outfile,live_u,live_v,live_logl)


# filename='updated_nlive_39.npz'
# live_points=numpy.load(filename)
# live_u=live_points['arr_0']
# live_v=live_points['arr_1']
# live_logl=live_points['arr_2']

# print(live_v[499])
# print(live_logl[499])