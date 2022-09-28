from os import listdir
from os.path import isfile,join

ganymede_files=[f for f in listdir('output_files/ganymede/')]
stampede_files=[f for f in listdir('output_files/stampede/')]

rm_ganymede=['MCMCA.out','MCMCU.out','temp.out']
rm_stampde=['slurm_temp.out','temp.out']

g=[item for item in ganymede_files if item not in rm_ganymede]
s=[item for item in stampede_files if item not in rm_stampde]

def get_systems(l):

    S=[]
    for item in l:
        f=item.split('.')
        if f[0] not in S:S.append(f[0])
    return S
        
ganymede_systems=get_systems(g)
stampede_systems=get_systems(s)

missing_systems=[g for g in ganymede_systems if not g in stampede_systems]
print(missing_systems)

systems=[g for g in ganymede_systems if g in stampede_systems]
systems.sort(key=int)
print(len(systems),systems)
# s=['1', '8', '12', '13', '17', '20', '25', '28', '32', '36', '39', '43', '44', '47', '48', '50', '54', '56', '57', '67', '70', '73', '76', '79', '80', '81', '83', '84', '85', '86', '88', '92', '93', '94', '95', '96', '106', '109', '120', '123', '126', '137']

with open('combined_analysis.txt','w') as fc:
    fc.write('filename\tstart_iter\thighest_value\thighest_m\tnext_m\tnext_m\n ')
    for system in systems:
        for c in ['ganymede','stampede']:
            for k in range(1,6):
                chain_filename=f'{system}_{c}_{k}.txt'
                with open(chain_filename,'r') as fs:
                    for i,lines in enumerate(fs):
                        y=lines.split()
                        if i==0:
                            high_value=y[0]
                            high_number=y[1]
                        elif i==1:
                            second_high=y[1]
                        elif i==2:
                            third_high=y[1]
                        else:break

                with open(f'output_files/{c}/{system}.{k}.out','r') as f:
                    for i,lines in enumerate(f):
                        x=lines.split()
                        if i==2:
                            it=x[2]
                            break
                    fc.write(f'{chain_filename}\t{it}\t{high_value}\t{high_number}\t{second_high}\t{third_high}\n')
