import numpy

s=['1', '8', '12', '13', '17', '20', '25', '28', '32', '36', '39', '43', '44', '47', '48', '50', '54', '56', '67', '70', '73', '76', '79', '80', '81', '83', '84', '85', '86', '88', '92', '93', '94', '95', '96', '106', '109', '120', '123', '126', '137']

with open('SpinlogQCatalog_el0.4.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        teff=float(x[2])
        teff_error=float(x[3])
        teff_lower=teff-teff_error
        teff_upper=teff+teff_error
        if numpy.logical_and(teff>4800,teff<5400):
            if x[0] in s:
                print(lines)



# Meibom and Mathieu 2005 trend for e-p enevelope 