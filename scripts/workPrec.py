from numpy import *
from time import *
import numpy.random as rnd
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Ellipse

WP = loadtxt('out/workprec.txt', skiprows=1)

P = WP[:,1]
Pp1 = P[1:]
Pm1 = P[:-1]
rat = Pp1/Pm1
print(rat)


fig = plt.figure(0)
fig.set_size_inches(7, 5, forward=True)
plt.rcParams.update({'font.size': 12})
plt.rcParams['xtick.major.size'] = 4.5
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.size'] = 4.5
plt.rcParams['ytick.major.width'] = 2
plt.semilogy(WP[:,0],WP[:,1], linewidth=3, linestyle='-', color='green')
plt.title(r'FLOPS-precision for free energy', size=18)
plt.xlabel('total FLOPS', size=14)
plt.ylabel('relative error', size=14)
plt.savefig('FEflops.eps', format='eps', dpi=2000)

fig = plt.figure(1)
fig.set_size_inches(7, 5, forward=True)
plt.rcParams.update({'font.size': 12})
plt.rcParams['xtick.major.size'] = 4.5
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.size'] = 4.5
plt.rcParams['ytick.major.width'] = 2
plt.semilogy(WP[:,1], linewidth=3, linestyle='-', color='blue')
plt.title(r'Convergence for free energy', size=18)
plt.xlabel('order of expansion (maximum n)', size=14)
plt.ylabel('relative error', size=14)
plt.savefig('FEconv.eps', format='eps', dpi=2000)

plt.show()
