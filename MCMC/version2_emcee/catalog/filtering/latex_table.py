
import csv
from email import header

data = []
with open('nominal_value_catalog_Iconv_cutoff.txt', 'r') as f:
    next(f)
    for lines in f:
        x = lines.split()
        val_list = ["%.4f" % float(val) for val in x[2:]]
        data.append([x[1],] + val_list)


header = ['KIC', r"$P_orb$ (days)", r"$P_{\star}$ (days)", r"e", r"Z", r"\$tau$ (Gyr)", r"$M_{\text{p}}$ ($M_{\odot}$)", r"$M_{\text{s}}$ ($M_{\odot}$)"]

with open('catalog.csv', 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(header)
    writer.writerows(data) 