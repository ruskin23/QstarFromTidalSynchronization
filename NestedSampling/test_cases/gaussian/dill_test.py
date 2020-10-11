import sys
import dill
import numpy

import matplotlib.pyplot  as plt

import dynesty
from dynesty import plotting as dyplot


def summary_plots(custom_results,
                  complete_results,
                  internal_samples,
                  truths,labels):


    fig, axes = plt.subplots(7,7,figsize=(25,25))

    fg,ax=dyplot.cornerplot(internal_samples, color='red',
                               truths=truths,labels=labels,truth_color='black',
                               show_titles=True, title_kwargs={'y': 1.05},
                               quantiles=None, fig=(fig, axes))


    fg.tight_layout()

    plt.show()
    #plt.savefig('custom.pdf')

with open('initial_samples_end.dill','rb') as f:
    break_samples=dill.load(f)
    break_samples=break_samples.results
with open('initial_samples.dill','rb') as f:
    complete_samples=dill.load(f)
    complete_samples=complete_samples.results
with open('internal_sampler.dill','rb') as f:
    internal_samples=dill.load(f)
    internal_samples=internal_samples.results

system_number=sys.argv[1]
with open('SpinlogQCatalog_el0.4.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        at_system=x[0]
        if system_number==at_system:
            feh_value=float(x[4])
            Porb_value=float(x[6])
            eccentricity_value=float(x[8])
            Pspin_value=float(x[12])
            break

truths=numpy.array([Porb_value,feh_value,eccentricity_value,((2*numpy.pi/14+2*numpy.pi/1.4)/2),8.0,0.8,4.5])
print(truths)
labels=[r"$Porb$",r"$feh$",r"eccentricity",r"$Wdisk$",r"$logQ$",r"$mass$",r"$age$"]
summary_plots(break_samples,complete_samples,internal_samples,truths,labels)
