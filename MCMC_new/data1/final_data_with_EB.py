def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

array1=[]
array2=[]

f3 = open('possible_candidates_data.txt','w')

with open('possible_candidates_KIC_simbad.txt','r') as f1:
    for line1 in f1:
        x = line1.strip()
        name1=''
        for v in x:
            if isfloat(v):
                #print(v)
                name1=name1+v
        #print(name1)
        with open('extract_useful_data.txt','r') as f2:
            for line2 in f2:
                x = line2.strip().split('\t')
                name2 =  x[0]
                if name1==name2:
                    f3.write(line2)



f3.close()
