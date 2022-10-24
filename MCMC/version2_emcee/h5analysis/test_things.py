from turtle import numinput
import numpy as np
import matplotlib.pyplot as plt

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


state = np.random.get_state()
random_state = np.random.RandomState()
random_state.set_state(state)
print(random_state.get_state())
random_state.seed()
print(random_state.get_state())

