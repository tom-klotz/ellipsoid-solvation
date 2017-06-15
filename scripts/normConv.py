from numpy import *
from time import *
import numpy.random as rnd
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Ellipse

###############
### USED WITH OUTPUT FROM EX3
###############

WP = loadtxt('out/normWorkPrec.txt', skiprows=1)
WP2 = loadtxt('../outperm/normWorkPrecSE.txt', skiprows=1) ##FIX LATER
WP3 = loadtxt('../outperm/normWorkPrecERF.txt', skiprows=1) ##FIX LATER
INT1 = loadtxt('out/normInt1Prec.txt', skiprows=1)
INT2 = loadtxt('out/normInt2Prec.txt', skiprows=1)
INT3 = loadtxt('out/normInt3Prec.txt', skiprows=1)
INT4 = loadtxt('out/normInt4Prec.txt', skiprows=1)

#E1INT1 = loadtxt('out/e1normInt1Prec.txt', skiprows=1)
#E1INT2 = loadtxt('out/e1normInt2Prec.txt', skiprows=1)
#E1INT3 = loadtxt('out/e1normInt3Prec.txt', skiprows=1)
#E1INT4 = loadtxt('out/e1normInt4Prec.txt', skiprows=1)

#E2INT1 = loadtxt('out/e2normInt1Prec.txt', skiprows=1)
#E2INT2 = loadtxt('out/e2normInt2Prec.txt', skiprows=1)
#E2INT3 = loadtxt('out/e2normInt3Prec.txt', skiprows=1)
#E2INT4 = loadtxt('out/e2normInt4Prec.txt', skiprows=1)

#E3INT1 = loadtxt('out/e3normInt1Prec.txt', skiprows=1)
#E3INT2 = loadtxt('out/e3normInt2Prec.txt', skiprows=1)
#E3INT3 = loadtxt('out/e3normInt3Prec.txt', skiprows=1)
#E3INT4 = loadtxt('out/e3normInt4Prec.txt', skiprows=1)


fig1 = plt.figure(0)
fig1.set_size_inches(7,5, forward=True)
plt.rcParams.update({'font.size': 12})
plt.rcParams['xtick.major.size'] = 4.5
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.size'] = 4.5
plt.rcParams['ytick.major.width'] = 2
plt.semilogy(WP[:,1],WP[:,2], linewidth=3, linestyle='-', color='green')
plt.title(r'FLOP-precision for $\gamma_{21}^{22}$', size=18)
plt.xlabel('total FLOPS', size=14)
plt.ylabel('relative error', size=14)
plt.savefig('figs/flopprec.eps', format='eps', dpi=2000)


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
plt.title(r'Convergence for $\gamma_{7}^{4}$', size=18)
plt.xlabel('quadrature nodes per integral', size=14)
plt.ylabel('relative error', size=14)
plt.savefig('figs/pointpreccomp.eps', format='eps', dpi=2000)

fig3 = plt.figure(2)
fig3.set_size_inches(7,5, forward=True)
plt.rcParams.update({'font.size': 12})
plt.rcParams['xtick.major.size'] = 4.5
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.size'] = 4.5
plt.rcParams['ytick.major.width'] = 2
plt.semilogy(INT1[:,0], INT1[:,1], linewidth=1.5, linestyle='-', marker='s', markersize=8, markevery=3, markerfacecolor='none', markeredgewidth=1, color='black')
plt.semilogy(INT2[:,0], INT2[:,1], linewidth=1.5, linestyle='-', marker='^', markersize=8, markevery=2, markerfacecolor='none', markeredgewidth=1, color='black')
plt.semilogy(INT3[:,0], INT3[:,1], linewidth=1.5, linestyle=':', color='black')#marker='d', markersize=8, markevery=3, markerfacecolor='none', markeredgewidth=1, color='black')
plt.semilogy(INT4[:,0], INT4[:,1], linewidth=1.5, linestyle='--', color='black')#marker='x', markersize=8, markevery=3, markerfacecolor='none', markeredgewidth=1, color='black')
plt.title(r'Tanh-sinh convergence for $\gamma_{7}^{4}$ (each integral)', size=18)
plt.xlabel('tanh-sinh nodes per integrand', size=14)
plt.ylabel('relative error', size=14)
plt.legend(['I_4', 'I_2', 'I_1', 'I_3'])
plt.savefig('figs/intConvComp.eps', format='eps', dpi=2000)

fig4 = plt.figure(3)
fig4.set_size_inches(7,5, forward=True)
plt.rcParams.update({'font.size': 12})
plt.rcParams['xtick.major.size'] = 4.5
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.size'] = 4.5
plt.rcParams['ytick.major.width'] = 2
plt.semilogy(WP[:,0],WP[:,2], linewidth=3, linestyle='-', color='black')
plt.title(r'Tanh-sinh convergence for $\gamma_{7}^{4}$', size=18)
plt.xlabel('tanh-sinh nodes per integral', size=14)
plt.ylabel('relative error', size=14)
plt.savefig('figs/pointprec.eps', format='eps', dpi=2000)

#fig5 = plt.figure(4)
#fig5.set_size_inches(7,5, forward=True)
#plt.rcParams.update({'font.size':12})
#plt.rcParams['xtick.major.size'] = 4.5
#plt.rcParams['xtick.major.width'] = 2
#plt.rcParams['ytick.major.size'] = 4.5
#plt.rcParams['ytick.major.width'] = 2
#plt.semilogy(E1INT4[:,0], E1INT4[:,1], linewidth=3, linestyle='-', color='blue')
#plt.semilogy(E2INT4[:,0], E2INT4[:,1], linewidth=3, linestyle='-', color='green')
#plt.semilogy(E3INT4[:,0], E3INT4[:,1], linewidth=3, linestyle='-', color='black')
#plt.xlabel('tanh-sinh nodes per integral', size=14)
#plt.ylabel('relative error', size=14)
#plt.savefig('pointprece1e2e3.eps', format='eps', dpi=2000)
plt.show()
