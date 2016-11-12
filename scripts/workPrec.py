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
#print(P)
#print(Pp1)
#print(Pm1)



plt.semilogy(WP[:,0],WP[:,1])


plt.show()
