from scipy import stats
from pathlib import Path
from directories import directories


system_number='10031409'
spin=4.039649675385522

home_dir=str(Path.home())
path=directories(home_dir)


observed_spin=dict()
with open(path.current_directory+'/catalog/filtering/Lurie_binaries_with_p1.txt','r') as f:
    for lines in f:
        x=lines.split()
        KIC=x[0]
        if KIC==system_number:
            observed_spin['value']=float(x[5])
            observed_spin['sigma']=abs(float(x[5])-float(x[6]))
            break
L=stats.norm(observed_spin['value'],observed_spin['sigma']).logpdf(spin)
print(L)