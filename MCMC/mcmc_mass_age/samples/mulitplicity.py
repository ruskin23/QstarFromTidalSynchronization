
#numbers =['113', '134', '133', '96', '33', '117', '112', '36', '79', '13', '57', '20', '2', '84', '32', '80', '30', '128', '123', '110', '7', '82', '119', '111', '121', '55', '65', '91', '14', '103', '130', '50', '25', '99', '137', '1', '86', '37', '132', '126', '139', '69', '51', '70', '76', '88', '115', '63', '142', '75', '105', '12', '67', '43', '108', '77', '38', '120', '28', '44', '83', '54', '68', '107', '85', '61', '94', '41', '17', '26', '29', '93', '131', '62', '90', '8', '109', '114', '4', '27', '95', '106', '35', '15', '92', '118', '49', '40', '31', '48', '73', '124', '72', '129', '81', '125', '23', '141', '59', '10', '52', '101', '78', '140', '56', '39', '18', '138', '100', '104', '5', '47']


numbers=['10']
for n in numbers:
    print(n)
    sample_file_name='MassAgeFehSamples_'+n+'.txt'
    with open('updated_samples/test_'+sample_file_name,'w') as f:
        f.write('mass'+'\t'+
                'age'+'\t'+
                'feh'+'\t'+
                'mulitplicity'+'\n')

    multiplicity=0
    with open('../'+sample_file_name,'r') as f:
        for i,lines in enumerate(f):
            x=lines.split()
            if i==0:previous_mass=float(x[1])
            current_mass=float(x[1])
            if current_mass==previous_mass:
                multiplicity=multiplicity+1
                paramter_set=x
            else:
                with open('updated_samples/test_'+sample_file_name,'a') as fw:
                    fw.write(paramter_set[1]+'\t'+
                            paramter_set[2]+'\t'+
                            paramter_set[3]+'\t'+
                            str(multiplicity)+'\n')
                previous_mass=float(x[1])
                multiplicity=1
