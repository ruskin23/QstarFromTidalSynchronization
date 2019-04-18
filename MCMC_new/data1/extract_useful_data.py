#Extract KIC and data with period less than 10, temperatures more than 5500 and
#less than 6200

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

def write_output(array1,array2):

    with open('extract_useful_data.txt','a') as f:
        for index,value in enumerate(array1):
            if index<7:f.write(value+'\t')
        f.write(array2[1] + '\n')

def create_array(lines):
    array=[]
    xd = lines.strip()
    check=''
    for value in xd:
        if value==' ':
            if check!='':array.append(check)
            check=''
            continue
        else: check=check+value
    return array

fname1 = 'eccentricity_teff_data.txt'
fname2='period_data.txt'

with open('extract_useful_data.txt','w') as f:
    f.write('#KIC' + '\t' + 'Teff1' + '\t' + 'Teff2' + '\t' + 'q' + '\t' + 'r1' +
            '\t' + 'r2' + '\t' + 'e' + '\t' + 'Per' + '\n')

f3=open('KIC_for_query_simbad.txt','w')
with open(fname1,'r') as f1, open(fname2,'r') as f2:
    for line1,line2 in zip(f1,f2):
        array1=create_array(line1)
        array2=create_array(line2)
        if  ('B' not in array1 and
            'B' not in array2 and
            float(array2[1])<10.0 and
            float(array1[1])>4300 and
            float(array1[1])<6200 and
            float(array1[2])>4300 and
            float(array1[2])<6200):

            print(float(array2[1]))
            write_output(array1,array2)
            f3.write('KIC'+array1[0]+'\n')



f3.close()
