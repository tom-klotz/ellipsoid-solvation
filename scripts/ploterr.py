from numpy import *
from time import *
import numpy.random as rnd
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Ellipse

X = loadtxt('out/xVals.txt')
Y = loadtxt('out/yVals.txt')
Xg,Yg = meshgrid(X, Y)
Z0 = transpose(loadtxt('out/solErr0.txt'))
Z1 = transpose(loadtxt('out/solErr1.txt'))
Z2 = transpose(loadtxt('out/solErr2.txt'))
Z3 = transpose(loadtxt('out/solErr3.txt'))
Z4 = transpose(loadtxt('out/solErr4.txt'))
Z5 = transpose(loadtxt('out/solErr5.txt'))
Z6 = transpose(loadtxt('out/solErr6.txt'))
Z7 = transpose(loadtxt('out/solErr7.txt'))

abc = loadtxt('out/otherinfo.txt')
chg = loadtxt('out/chargeXYZ.txt')
col = loadtxt('out/chargeMag.txt')
ellAlpha = .3
color = rnd.rand(3)
ells = [Ellipse(xy=[0,0], width=2*abc[0], height=2*abc[1], angle=0) for i in range(8)]
area = 50


fig0 = plt.figure(0)
ax0 = fig0.gca(projection='3d')
ax0.plot_surface(Xg, Yg, Z5)


fig = plt.figure(1)

ax = fig.add_subplot(221, aspect='equal')
ax.add_artist(ells[5])
#ax.plot_surface(Xg,Yg,Z5)
ax.scatter(chg[:,0],chg[:,1], s=area, c=col, alpha=0.5)
ells[5].set_alpha(ellAlpha)
ells[5].set_facecolor(color)
ax.set_title('error 5')


plt.show()
