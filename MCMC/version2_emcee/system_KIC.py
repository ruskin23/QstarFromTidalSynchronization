
d=[]
with open('/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/version2_emcee/catalog/filtering/nominal_value_catalog_temp_cutoff.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        d.append(x[1])


print(len(d))
print(*d, sep='\t')

#10264202        10330495        10385682        10454401        10935310        11147276        11228612        11232745        11233911        11252617        11403216        11616200        11704044        12470530        2305543 3241344     3348093 3973504 4352168 4678171 5730394 7732791 8356054 8543278 8580438 8618226 8957954 10091257        10257903        10711913        10879213        10936427        10960995        10992733        11200773        11391181   12004679 12351927        2447893 2449090 2852560 2860788 3098194 3248019 3323289 3344427 3834364 3838496 4252226 4285087 4346875 4380283 4757331 4773155 4902030 5022440 5300878 5393558 5652260 5802470 5871918 6029130 6128027 6283224 6312521     6359798 6449552 6863840 6927629 6949550 6962018 7025540 7362852