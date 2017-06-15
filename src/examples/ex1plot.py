from numpy import *
from time import *
import numpy.random as rnd
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Ellipse


# First column should be ellipsoidal inside (converging)
# Second column should be spherical inside (diverging)
# Third column should be ellipsoidal outside (converging)
# Fourth column should be spherical outside (converging slower)
D = loadtxt('out/ex1.txt')
x = range(1,26)

fig0 = plt.figure(0)
fig0.set_size_inches(5, 2.5, forward=True)
plt.rcParams.update({"font.size": 8})
plt.semilogy(x, D[:,0], linestyle='-', color='black', linewidth=2)
plt.semilogy(x, D[:,1], linestyle='--', color='black', linewidth=2)
plt.legend(['Ellipsoidal expansion', 'Spherical expansion'])
plt.title('Coulomb Convergence Inside Reference Ellipsoid')
plt.xlabel('Maximum expansion order (N)')
plt.ylabel('Relative error')
plt.tight_layout()
plt.savefig('figs/CoulCompInside.eps', format='eps', dpi=2000)


fig1 = plt.figure(1)
fig1.set_size_inches(5, 2.5, forward=True)
plt.rcParams.update({"font.size": 8})
plt.semilogy(x, D[:,2], linestyle='-', color='black', linewidth=2)
plt.semilogy(x, D[:,3], linestyle='--', color='black', linewidth=2)
plt.legend(['Ellipsoidal expansion', 'Spherical expansion'])
plt.xlabel('Maximum expansion order (N)')
plt.title('Coulomb Convergence Outside Reference Ellipsoid')
plt.ylabel('Relative error')
plt.tight_layout()
plt.savefig('figs/CoulCompOutside.eps', format='eps', dpi=2000)
plt.show()
