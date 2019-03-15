import pickle

#data = [1,2,3,4]
data = 1
with open('data_test.pcikle','wb') as f :
    pickle.dump(data,f)



with open('data_test.pcikle','rb') as f:
    x = pickle.load(f)
    print(x)
