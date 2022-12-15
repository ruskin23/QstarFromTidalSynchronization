import numpy as np
import matplotlib.pyplot as plt
import pickle


# Fixing random state for reproducibility
np.random.seed(19680801)

# some random data
x = np.random.randn(1000)
y = np.random.randn(1000)


def scatter_hist(x, y, ax, ax_histy):
    # no labels
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(x, y)
    ax.set_ylabel(r'$\log_{10}{Q_{\ast}^{\prime}}$')
    ax.set_xlabel(r'$\log_{10}{P_{tide}}$')

    # now determine nice limits by hand:
    binwidth = 0.25
    xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
    lim = (int(xymax/binwidth) + 1) * binwidth

    bins = np.arange(-lim, lim + binwidth, binwidth)
    ax_histy.hist(y, bins=bins, orientation='horizontal')

    
    plt.savefig('test.png')




def test_dict(dict_a):

    dict_a['2'] = dict_a['2'][5:]
# Start with a square Figure.
# fig = plt.figure(figsize=(6, 6))
# Add a gridspec with two rows and two columns and a ratio of 1 to 4 between
# the size of the marginal axes and the main axes in both directions.
# Also adjust the subplot parameters for a square plot.
# gs = fig.add_gridspec(1, 2,  width_ratios=(4, 1), wspace=0.05)
# Create the Axes.
# ax = fig.add_subplot(gs[0, 0])
# ax_histy = fig.add_subplot(gs[0, 1], sharey=ax)
# Draw the scatter plot and marginals.
# scatter_hist(x, y, ax, ax_histy)


# state = np.random.get_state()
# random_state = np.random.RandomState()
# random_state.set_state(state)
# print(random_state.get_state())
# random_state.seed()
# print(random_state.get_state())



# a = dict()
# a['1'] = 4
# a['2'] = np.array([1,2,3,4,5,6,7,8,9,10])

# print(a)
# test_dict(a)
# print(a)


# a = np.random.rand(4,4)
# print(a)
# b = np.random.rand(4,4) - 0.5
# print(b)

# c = np.empty(np.shape(a))

# c[b<0] = a[b<0]*10
# c[b>0] = a[b>0]*0.1

# print(c)

# print(b<0)

# a = np.random.rand(4,4)
# print(a)
# for i,b in enumerate(a):
#     for j in range(len(b)):
#         print(a[i,j])


# test_tup = [(1,2), (3,4), (5,6), (7,8), (9,10)]
# f = open('test.pickle', 'wb')
# pickle.dump(test_tup, f)
# f.close()


# f = open('test.pickle', 'rb')
# test_tup = pickle.load(f)
# print(test_tup)
# f.close()

# f = open('sampled_params.pickle', 'rb')
# tup_array = pickle.load(f)
# print(tup_array)


# x = np.arange(4)
# y = 4

# xx, yy = np.meshgrid(x, y)
# print(xx)
# print(yy)

# diff = np.squeeze(xx - yy )
# print(diff)

# avg = np.average(diff, axis = -1)
# print(avg)

# a = np.ones((2,4,4), dtype = float)
# b = np.arange(16).reshape(4,4)
# a[0] = b
# print(a[1])

# def equal_vowels(arr):

#     vowels = ['a', 'e', 'i', 'o', 'u', 'A', 'E', 'I', 'O', 'U']

#     upper = arr[:int(len(arr)/2)]
#     upper_count = 0

#     lower = arr[int(len(arr)/2):]
#     lower_count = 0

#     for i in range(int(len(arr)/2)):

#         if upper[i] in vowels: upper_count += 1
#         if lower[i] in vowels: lower_count += 1
    
#     return upper_count == lower_count

# def call_equal_vowel():
#     import sys
#     import random
#     import string

#     len_random = int(sys.argv[1])

#     word = ''.join(random.choices(string.ascii_uppercase + string.digits, k=len_random))
#     print(word)
#     # word =sys.argv[1]
#     print(equal_vowels(word))