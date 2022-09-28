import numpy


def get_parameters_logfile(logfilename):

    _simple_quantities=['primary_mass',
                        'secondary_mass',
                        'feh',
                        'age',
                        'Wdisk',
                        'orbital_period',
                        'eccentricity',
                        'phase_lag_max']

    parameters=dict()

    with open('files/'+logfilename,'r') as f:
        for i,lines in enumerate(f):

            if i>6 and i<42 and i%2==1:

                x=lines.split()

                if x[0] in _simple_quantities:
                    parameters[x[0]]=float(x[1])

                if x[0]=='tidal_frequency_breaks':
                    x=lines.split()
                    print(len(x))
                    if len(x)==2:
                        parameters[x[0]]=numpy.atleast_1d(float(x[1][1:-1]))
                    elif len(x)==3:
                        a=float(x[1][1:])
                        b=float(x[-1][:-1])
                        parameters[x[0]]=numpy.array([a,b])
                    elif len(x)==4:
                        a=float(x[1][1:])
                        b=float(x[-2])
                        parameters[x[0]]=numpy.array([a,b])

                if x[0]=='tidal_frequency_powers':
                    if len(x)==4:
                        a=float(x[2])
                        b=float(x[3][0:-1])
                        value=numpy.array([a,b])
                    elif len(x)==5:
                        a=float(x[1][1:])
                        b=float(x[2])
                        c=float(x[3])
                        value=numpy.array([a,b,c])
                    parameters[x[0]]=numpy.array(value)
    return parameters

p=get_parameters_logfile('10215422_20220824151945_575947.log')
print(p)
tp=p['tidal_frequency_breaks']
print(len(tp))