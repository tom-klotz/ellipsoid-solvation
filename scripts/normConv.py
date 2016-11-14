from numpy import *
from time import *
import numpy.random as rnd
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Ellipse

WP = loadtxt('out/normWorkPrec.txt', skiprows=1)

plt.figure(0)
plt.semilogy(WP[:,1],WP[:,2])
plt.figure(1)
plt.semilogy(WP[:,0],WP[:,2])
plt.show()
