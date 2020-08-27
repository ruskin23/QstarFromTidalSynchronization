


old=[]
new=[]

with open('SpinlogQCatalog_el0.4.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        old.append(x[0])


with open('WindeCatalog.txt','r') as f:
    next(f)
    for lines in f:
        x=lines.split()
        new.append(x[0])



running=['85', '76', '96', '81', '80', '36', '83', '84', '94', '32', '106',
         '123', '50', '39',
         '56', '126', '54', '70', '88', '67', '95', '25', '137', '1', '86',
         '43', '73', '92',
         '93', '79', '47', '109', '44', '48', '17', '8', '12', '20', '57',
         '120', '28',
         '13']

print('Common between running and new:')
print('number systems',len(set(new) & set(running)))
print('system numbers:',set(new) & set(running))
print('remaining systems:',set(new)-set(running))

print(len(set(old) & set(new)))
print(set(old) & set(new))


