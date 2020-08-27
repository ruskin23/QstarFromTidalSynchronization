import numpy as np
import pandas as pd
import sys
from sklearn import preprocessing
from sklearn.multioutput import MultiOutputRegressor
from sklearn.ensemble import RandomForestRegressor,AdaBoostRegressor,GradientBoostingRegressor
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import train_test_split

import matplotlib.pyplot as plt
import pickle

print("Read Data")

df = pd.read_csv("extracted.csv",names = ['Porb','ecc','Wdisk','logQ','Mass_P','age','feh','Mass_S','Porb_i','ecc_i'])

'''
labels are Porb_i and ecc_i
'''
#print(df)

Y = df[['Porb_i', 'ecc_i']].copy()

X = df.drop('ecc_i',axis = 1)
X = X.drop('Porb_i',axis = 1)

min_max_scaler = preprocessing.MinMaxScaler()
x_scaled = min_max_scaler.fit_transform(X)
X_norm = pd.DataFrame(x_scaled)
#"Normalized X"

#print(X_norm)

X_trn, X_tst, y_trn, y_tst = train_test_split(X_norm, Y,
                                                  random_state=4,
                                                  test_size=0.1, shuffle=True)

#print(len(X_trn),len(X_tst))

X_train, X_val, y_train, y_val = train_test_split(X_trn, y_trn,
                                                  random_state=4,
                                                  test_size=0.2, shuffle=True)
print("Training Samples: ",len(X_train),"Validation Set: ",len(X_val), "Test set: ", len(X_tst))

losses = []
#for depth in range(2,20):

regr = RandomForestRegressor(n_estimators=100, max_depth=7,
                             random_state=0, criterion='mae' ,n_jobs=10)

#regr = AdaBoostRegressor(n_estimators=500, random_state=0, loss='square' )
clf = MultiOutputRegressor(regr).fit(X_train,y_train)

# scores = cross_val_score(clf, X_trn, y_trn, cv=5, scoring='neg_root_mean_squared_error')
# print(scores)
predictions = clf.predict(X_val)
print(predictions)
#f = open("predictions.txt",mode='w+')
y_out = np.array(y_val)
# for i in range(len(predictions)):
#     #temp = [predictions[i] , y_tst[i]]
#     f.write(repr(predictions[i] - y_tst[i])+'\n')
loss = np.sum((predictions - y_out)**2)
print("For the depth = ", 16, " loss is: ",loss)
losses.append(loss)
pred_test = clf.predict(X_tst)
lossontest = np.sum((pred_test - np.array(y_tst))**2)
print("for the test set the loss is: ",lossontest)

with open('model_16.pickle','wb') as f:
    pickle.dump(clf,f)


#plt.figure()
#plt.title("Depthwise analysis")
#plt.xlabel(" Tree Depths ")
#plt.ylabel("Squared Loss")
#plt.semilogy(np.arange(2,20),np.array(losses))
#plt.show()
