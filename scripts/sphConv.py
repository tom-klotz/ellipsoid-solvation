from numpy import *
from time import *
import numpy.random as rnd
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Ellipse

WP = loadtxt('out/sphConv.txt', skiprows=1)

fig = plt.figure(0)
fig.set_size_inches(7, 5, forward=True)
plt.rcParams.update({'font.size': 12})
plt.rcParams['xtick.major.size'] = 4.5
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.size'] = 4.5
plt.rcParams['ytick.major.width'] = 2
plt.loglog(WP[:,0], WP[:,1], linewidth=3, marker='o', mew=2, ms=10)
plt.title('Converge towards Born ion', size=18)
plt.xlabel(r'Axis perturbation size ($\Delta$)', size=14)
plt.ylabel('Deviation from Born ion (relative error)', size=14)
plt.savefig('sphConv.eps', format='eps', dpi=2000)
plt.show()
