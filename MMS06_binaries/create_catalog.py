



with open('final_catalog.txt','w') as f:
    f.write('system'+'\t'+
            'Porb'+'\t'+
            'e'+'\t'+
            'sigma_e'+'\t'+
            'primary_mass'+'\t'+
            'mass_function'+'\t'+
            'P_star'+'\t'+
            'Pstar_error'+'\n')


    with open('Table1.txt','r') as f1:
        next(f1)
        for lines in f1:
            x=lines.split()
            system=x[0]
            Porb=float(x[4])
            e=float(x[5])
            sigma_e=float(x[6])
            primary_mass=float(x[9])
            P_star=float(x[7])
            Pstar_error=float(x[8])
            with open('rv_fit.txt','r') as f2:
                next(f2)
                for lines in f2:
                    y=lines.split()
                    if y[0]==system:
                        K=float(y[1])
                        ecc=float(y[2])
                        f_M=Porb*K*K*K*((1-ecc*ecc)**1.5)*1.03e-7
                        break
            f.write(system+'\t'+
                    repr(Porb)+'\t'+
                    repr(e)+'\t'+
                    repr(sigma_e)+'\t'+
                    repr(primary_mass)+'\t'+
                    repr(f_M)+'\t'+
                    repr(P_star)+'\t'+
                    repr(Pstar_error)+'\n')
