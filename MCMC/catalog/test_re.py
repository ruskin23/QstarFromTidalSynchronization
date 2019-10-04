import re

pattern = '9410241'
text = 'randomly wrting something 110'
with open('Meiborn_2011_NGC6811_rot_per.txt','r') as f:
    filetext=f.read()
    if (re.search(pattern,filetext)) : print(filetext)
