from imp import NullImporter
from re import T
import sys
import os
import random
import pandas
import numpy

from pathlib import Path
from directories import directories

home_dir=str(Path.home())
path=directories(home_dir)
sys.path.append(path.poet_path+'/PythonPackage')
sys.path.append(path.poet_path+'/scripts')


from stellar_evolution.library_interface import library

from stellar_evolution.manager import StellarEvolutionManager
from orbital_evolution.evolve_interface import library as \
    orbital_evolution_library
from stellar_evolution.derived_stellar_quantities import\
TeffK,\
LogGCGS,\
RhoCGS

serialized_dir = path.poet_path +  "/stellar_evolution_interpolators"
manager = StellarEvolutionManager(serialized_dir)
interpolator = manager.get_interpolator_by_name('default')



def stellar_parameter_check(interpolator,feh,t,m1,m2):

    
    if numpy.logical_and(numpy.logical_and(feh>-1.014,feh<0.537)
                        ,
                        numpy.logical_and(numpy.logical_and(m1>0.4,m1<1.2),
                                            numpy.logical_and(m2>0.4,m2<1.2))
                    ):
                    age_max_m1=interpolator('radius', m1, feh).max_age
                    return numpy.logical_and(t>8,(10**t)/1e9<age_max_m1)
    else: return False
                    # if numpy.logical_and(t>8,(10**t)/1e9<age_max_m1):
                    #     return True



if __name__=='__main__':
    spin_dict=dict()
    spin_KIC=[]
    with open('Lurie_binaries_with_p1.txt','r') as f:
        for lines in f:
            x=lines.split()
            spin_KIC.append(int(x[0]))
            spin_dict[x[0]]=dict(value=float(x[5]),sigma=abs(float(x[6])-float(x[5])))
    n = 0
    with open('windemuth_orbital_raw.txt', 'r') as f:
        next(f)
        next(f)
        for lines in f:
            x = lines.split()
            print(x[0])
            if x[0] in spin_KIC:n+=1
    print(n)

    # win_KIC=[]
    # with open('windemuth_orbital_raw.txt','r') as f:
    #     next(f)
    #     next(f)
    #     for lines in f:
    #         x=lines.split()
    #         win_KIC.append(int(x[0]))

    # print(set(spin_KIC) - set(win_KIC))

    # serialized_dir =  "/home/ruskin/projects/poet/stellar_evolution_interpolators"
    # manager = StellarEvolutionManager(serialized_dir)
    # interpolator = manager.get_interpolator_by_name('default')

    # eccentricity_path=os.path.join('/home/ruskin/projects/poet','eccentricity_expansion_coef.txt').encode('ascii')

    # orbital_evolution_library.read_eccentricity_expansion_coefficients(
    # eccentricity_path)

    KIC_not_in_win=[]
    fo=open('filtered_spin_catalog.txt','w')
    fo.write('Number'+'\t'+
            'KIC'+'\t'+
            'Porb'+'\t'+
            'esinw'+'\t'+
            'esinw+sigma'+'\t'+
            'esinw-sigma'+'\t'+
            'ecosw'+'\t'+
            'ecos+sigma'+'\t'+
            'ecos-sigma'+'\t'+
            'z'+'\t'+
            'z+sigma'+'\t'+
            'z-sigma'+'\t'+
            'tau (log10 yr)'+'\t'+
            'tau+sigma'+'\t'+
            'tau-sigma'+'\t'+
            'm1'+'\t'+
            'm1+sigma'+'\t'+
            'm1-sigma'+'\t'+
            'm2'+'\t'+
            'm2+sigma'+'\t'+
            'm2-sigma'+'\t'+
            'spin'+'\t'+
            'sigma_spin'+'\n')
    n=0
    with open('windemuth_orbital_raw.txt','r') as fp:
        next(fp)
        next(fp)
        for i,lines1 in enumerate(fp):
            params=[]
            x=lines1.split()
            wind_KIC=int(x[0])
            if wind_KIC in spin_KIC:
                n=n+1
                porb=float(x[1])
                esinw=float(x[7])
                ecosw=float(x[10])
                eccentricity=numpy.sqrt(esinw**2 + ecosw**2)
                if numpy.logical_and(porb<60,porb>0) and eccentricity<0.50:
                    params.append(str(i+1))
                    params.append(x[0])
                    params.append(x[1])
                    params.append(x[7])
                    params.append(x[8])
                    params.append(x[9])
                    params.append(x[10])
                    params.append(x[11])
                    params.append(x[12])
                    with open('windemuth_stellar_raw.txt','r') as fs:
                        next(fs)
                        next(fs)
                        for lines2 in fs:
                            y=lines2.split()
                            if y[0]==x[0]:
                                z=float(y[1])
                                feh=library.feh_from_z(z)
                                t=float(y[4])
                                m1=float(y[7])
                                m2=float(y[10])
                                if stellar_parameter_check(interpolator,feh,t,m1,m2):
                                    params.append(y[1])
                                    params.append(y[2])
                                    params.append(y[3])
                                    params.append(y[4])
                                    params.append(y[5])
                                    params.append(y[6])
                                    params.append(y[7])
                                    params.append(y[8])
                                    params.append(y[9])
                                    params.append(y[10])
                                    params.append(y[11])
                                    params.append(y[12])
                                    params.append(str(spin_dict[x[0]]['value']))
                                    params.append(str(spin_dict[x[0]]['sigma']))
                                    line='\t'.join(params)
                                    fo.write(line)
                                    fo.write('\n')
    fo.close()


    f_nominal=open('nominal_value_catalog.txt','w')
    f_nominal.write('Number'+'\t'+
            'KIC'+'\t'+
            'Porb'+'\t'+
            'spin'+'\t'+
            'spin_error'+'\t'+
            'eccentricity'+'\t'+
            'feh'+'\t'+
            'age (Gyr)'+'\t'+
            'm1'+'\t'+
            'm2'+'\n')

    with open('filtered_spin_catalog.txt','r') as f_filtered_spin_catalog:
        next(f_filtered_spin_catalog)
        PO=[]
        PS=[]
        for lines in f_filtered_spin_catalog:
            x=lines.split()
            number=x[0]
            KIC=x[1]
            porb = float(x[2])
            # porb=round(float(x[2]), 3)
            PO.append(porb)
            esinw=float(x[3])
            ecosw=float(x[6])
            eccentricity=numpy.sqrt(esinw**2 + ecosw**2)
            feh=library.feh_from_z(float(x[9]))
            age=(10**(float(x[12])))/1e9
            m1=x[15]
            m2=x[18]
            spin = float(x[21])
            spin_error = float(x[22])
            # spin = round(float(x[21]), 3)
            # spin_error = round(float(x[22]), 3)
            PS.append(spin)
            f_nominal.write(f'{number}\t{KIC}\t{porb}\t{spin}\t{spin_error}\t{eccentricity}\t{feh}\t{age}\t{m1}\t{m2}\n')
    f_nominal.close()

    # Total Systems with morph<0.5
    maxlike_data = pandas.read_csv(
        'paper_maxlike_pars.dat',
        sep=r'\s+'
    )
    KIC = maxlike_data[maxlike_data['morph']<0.5]['#KIC'].to_numpy()

    #second nominal file where min(I(conv))>0
    with open('nominal_value_catalog_Iconv_cutoff.txt','w') as fnew:
        fnew.write('Number\tKIC\tPorb\tspin\teccentricity\tfeh\tage(Gyr)\tm1\tm2\n')
        with open('nominal_value_catalog.txt','r') as f:
            next(f)
            for lines in f:
                x=lines.split()
                if int(x[1]) in KIC:
                    feh=float(x[6])
                    age=float(x[7])
                    m1=float(x[8])
                    m2=float(x[9])
                    quantity_radius=interpolator('radius',m1, feh)
                    quantity_lum=interpolator('lum',m1, feh)
                    t_age=numpy.linspace(5e-3,age,1000)
                    # try:
                    T=TeffK(quantity_radius,quantity_lum)(age)
                    I=interpolator('iconv',m1,feh)
                    if min(I(t_age))>0:
                        # print(x[1])
                        fnew.write(lines)
                        # print(f'System {x[0]} Temp = {T} M1 = {m1} feh = {feh} age = {age}')
                    # except:
                    #     continue




# PO,PS=sorted(PO),sorted(PS)
# print(max(PO))
# for p,s in zip(PO,PS):
#     print(f'{p}\t{s}')

