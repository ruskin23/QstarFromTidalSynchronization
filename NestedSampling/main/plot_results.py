import sys
import dill
import numpy

import matplotlib.pyplot  as plt

import dynesty
from dynesty import plotting as dyplot


def summary_plots(partial_results,labels):

    fig, axes = plt.subplots(7,7,figsize=(25,25))


    fg,ax=dyplot.cornerplot(partial_results, color='red',
                               labels=labels,
                               show_titles=True, title_kwargs={'y': 1.05},
                               quantiles=None, fig=(fig, axes))


    fg1,ax1=dyplot.runplot(partial_results)
    fg.tight_layout()

    plt.show()
    #plt.savefig('custom.pdf')


system_number=sys.argv[1]
#sampler_file='/home/ruskin/projects/QstarFromTidalSynchronization/NestedSampling/main/results/initial_sampling_saved_'+system_number+'.dill'
sampler_file='/home/ruskin/projects/QstarFromTidalSynchronization/NestedSampling/main/results/kartof/initial_sampling_saved_76.dill'
with open(sampler_file,'rb') as f:
    partial_sampler=dill.load(f)
    partial_results=partial_sampler.results

samples=partial_results.samples
logl=partial_results.logl
print(logl.size)
print(max(logl))

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
summary_plots(partial_results,labels)
