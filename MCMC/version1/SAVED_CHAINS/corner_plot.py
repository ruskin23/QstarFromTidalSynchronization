import corner
import pickle
import numpy
import matplotlib.pyplot as plt
import sys


data=[]

system=sys.argv[1]

with open('complete_chains.pickle','rb') as f:
    D=pickle.load(f)

for system_name,parameters in D.items():
    if system_name==system:
        for param,value in parameters.items():
            data.append(value)

d=numpy.vstack(data)
d=d.T
# with open('AcceptedParameters.txt','r') as f:
#     next(f)
#     for lines in f:
#         x=lines.split()
#         data.append(float(x[1]))
#         data.append(float(x[2]))
#         data.append(float(x[3]))
#         data.append(float(x[4]))
#         data.append(float(x[5]))
#         data.append(float(x[6]))
#         data.append(float(x[7]))
#         data.append(float(x[15]))

# data=np.array(data)
# d=data.reshape([len(data)//8,8])

figure=corner.corner(d,
                     labels=[r"$Porb$",r"eccentricity",r"$Wdisk$",r"$logQ$",r"$mass$",r"$age$",r"$feh$",r"$Pspin$"],
                     color='k',
                     show_titles=True)



plt.figure(1)
plt.show()
# plt.savefig('CornerPlot'+system+'.eps')

# 48      1888    45.03737867335107       0.957310002880258       1.0121208649129803
# 70      942     16.610445445627587      0.5927594699570033      1.014241600861263
# 43      2061    95.19782011890291       1.2389218099604604      1.0182324420912452
# 88      2811    164.63670862366678      1.4712717419347696      1.0195354289776781
# 13      1546    17.979264198830155      0.27265102248442263     1.020787309086212
# 96      11470   685.2380114545765       1.2260616971000018      1.0240309421840805
# 44      1052    47.47930664563399       0.8623431983285093      1.0253713419209924
# 86      9342    457.52333783486654      0.8745685733832926      1.0275659880598216
# 123     1540    166.16227889140367      1.889429132864187       1.0278407301681118
# 17      1763    105.38069751301151      1.0359811575483646      1.0281684065832493
# 93      2128    102.51181291600886      0.8338291318659383      1.0282524572982437
# 94      1459    72.8021075721828        0.6623032885357267      1.0366560002954013
# 84      1003    183.80998937583772      1.6082669259107418      1.0549653393227998
# 67      9858    1310.5155586940252      1.1275942870280689      1.0572582219054585
# 50      5779    1129.281269315799       1.6143656230488848      1.0587124953057296
# 85      5868    610.0674023637724       0.6914555854971906      1.0724674847680624
# 1       12921   1973.3098452906474      0.9158843012706548      1.0801249016674521
# 28      806     340.7882399584458       2.2325510112212563      1.0900208411175178
# 95      16181   5730.572356076061       1.8247174726444217      1.0927147164205127
# 39      7343    337.86855138788735      0.21812229080563705     1.1003686759038023
# 36      2005    978.090096357701        1.964385921331959       1.117065842988158
# 8       1091    920.366337811549        3.2558297574213193      1.1216895786967307
# 25      589     226.73678127329953      1.325017940189067       1.1352656751000012
# 32      2435    647.309019719363        0.8943880254172116      1.1387779712227093
# 126     8986    5080.713992117027       1.8213130172920944      1.14469467747843
# 81      6232    1304.5188865017667      0.5886698454208001      1.1642297221661653
# 106     872     692.1523357068704       2.2158526889625927      1.1649329018129375
# 109     1753    1512.9182602394606      2.323527996968583       1.17084027785384
# 137     6346    1298.7232116514706      0.5224540435402247      1.179642247005865
# 47      397     194.97809699582425      1.0312209775184031      1.2139771543856055
# 76      2027    85.4915741194571        0.08700955549099143     1.2182937623491132
# 83      1336    563.8719347931662       0.8452394566844774      1.224168732171774
# 12      716     1083.4584394252352      2.647597201564769       1.2530539736319377
# 20      2246    1023.7064262305357      0.727111758295195       1.275306250935465
# 54      2721    505.41623150139753      0.29559524168521206     1.2759364086907434
# 92      2096    1339.4547470419475      0.8084035788521982      1.337921908641234
# 120     586     325.04404177541693      0.6617056713247935      1.355195717703587
# 56      4331    3009.099126764487       0.8004954327041272      1.3666413590918025
# 79      868     1200.7325684871796      1.1126928693276004      1.4973567690379235
# 73      2303    4160.160166326555       1.2748046047796622      1.5545334978760463
# 80      7051    11525.001661239556      0.39203551638217216     2.2735819127402728