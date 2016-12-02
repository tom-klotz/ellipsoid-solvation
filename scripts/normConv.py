from numpy import *
from time import *
import numpy.random as rnd
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Ellipse

WP = loadtxt('outperm/normWorkPrec.txt', skiprows=1)
WP2 = loadtxt('outperm/normWorkPrecSE.txt', skiprows=1)
WP3 = loadtxt('outperm/normWorkPrecERF.txt', skiprows=1)
INT1 = loadtxt('outperm/normInt1Prec.txt', skiprows=1)
INT2 = loadtxt('outperm/normInt2Prec.txt', skiprows=1)
INT3 = loadtxt('outperm/normInt3Prec.txt', skiprows=1)
INT4 = loadtxt('outperm/normInt4Prec.txt', skiprows=1)


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
plt.savefig('pointpreccomp.pdf', format='eps', dpi=2000)

fig3 = plt.figure(2)
fig3.set_size_inches(7,5, forward=True)
plt.rcParams.update({'font.size': 12})
plt.rcParams['xtick.major.size'] = 4.5
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.size'] = 4.5
plt.rcParams['ytick.major.width'] = 2
plt.semilogy(INT1[:,0], INT1[:,1], linewidth=3, linestyle='-', color='blue')
plt.semilogy(INT2[:,0], INT2[:,1], linewidth=3, linestyle='-', color='green')
plt.semilogy(INT3[:,0], INT3[:,1], linewidth=3, linestyle='-', color='red')
plt.semilogy(INT4[:,0], INT4[:,1], linewidth=3, linestyle='-', color='black')
plt.title(r'Tanh-sinh convergence for $\gamma_3^5$ (each integral)', size=18)
plt.xlabel('tanh-sinh nodes per integral', size=14)
plt.ylabel('relative error', size=14)
plt.savefig('intConvComp.pdf', format='eps', dpi=2000)

fig4 = plt.figure(3)
fig4.set_size_inches(7,5, forward=True)
plt.rcParams.update({'font.size': 12})
plt.rcParams['xtick.major.size'] = 4.5
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.size'] = 4.5
plt.rcParams['ytick.major.width'] = 2
plt.semilogy(WP[:,0],WP[:,2], linewidth=3, linestyle='-', color='blue')
plt.title(r'Tanh-sinh convergence for $\gamma_3^5$', size=18)
plt.xlabel('tanh-sinh nodes per integral', size=14)
plt.ylabel('relative error', size=14)
plt.savefig('pointprec.eps', format='eps', dpi=2000)

plt.show()
