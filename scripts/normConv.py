from numpy import *
from time import *
import numpy.random as rnd
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Ellipse

WP = loadtxt('out/normWorkPrec.txt', skiprows=1)
WP2 = loadtxt('out/normWorkPrecSE.txt', skiprows=1)
WP3 = loadtxt('out/normWorkPrecERF.txt', skiprows=1)

fig1 = plt.figure(0)
fig1.set_size_inches(7,5, forward=True)
plt.rcParams.update({'font.size': 12})
plt.rcParams['xtick.major.size'] = 4.5
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.size'] = 4.5
plt.rcParams['ytick.major.width'] = 2
plt.semilogy(WP[:,1],WP[:,2], linewidth=3, linestyle='-', color='green')
plt.title(r'FLOP-precision for $\gamma_3^5$', size=18)
plt.xlabel('total FLOPS', size=14)
plt.ylabel('relative error', size=14)
plt.savefig('flopprec.eps', format='eps', dpi=2000)


fig2 = plt.figure(1)
fig2.set_size_inches(7,5, forward=True)
plt.rcParams.update({'font.size': 12})
plt.rcParams['xtick.major.size'] = 4.5
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.size'] = 4.5
plt.rcParams['ytick.major.width'] = 2
plt.semilogy(WP[:,0],WP[:,2], linewidth=3, linestyle='-', color='blue')
plt.semilogy(WP3[:,0], WP3[:,2], linewidth=3, linestyle='-', color='green')
plt.semilogy(WP2[:,0], WP2[:,2], linewidth=3, linestyle='-', color='red')
plt.title(r'Convergence for $\gamma_3^5$', size=18)
plt.xlabel('quadrature nodes per integral', size=14)
plt.ylabel('relative error', size=14)
plt.savefig('pointpreccomp.eps', format='eps', dpi=2000)

fig3 = plt.figure(2)
fig3.set_size_inches(7,5, forward=True)
plt.rcParams.update({'font.size': 12})
plt.rcParams['xtick.major.size'] = 4.5
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.size'] = 4.5
plt.rcParams['ytick.major.width'] = 2
plt.semilogy(WP2[:,0],WP2[:,2], linewidth=3, linestyle='-', color='blue')
plt.title(r'ERF Convergence for $\gamma_3^5$', size=18)
plt.xlabel('# tanh-sinh nodes', size=14)
plt.ylabel('relative error', size=14)
plt.savefig('pointprec.eps', format='eps', dpi=2000)

plt.show()
