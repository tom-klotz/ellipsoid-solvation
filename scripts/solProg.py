from numpy import *
from time import *
import numpy.random as rnd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse


X = loadtxt('out/xVals.txt')
Y = loadtxt('out/yVals.txt')
Z0 = transpose(loadtxt('out/sol0.txt'))
Z1 = transpose(loadtxt('out/sol1.txt'))
Z2 = transpose(loadtxt('out/sol2.txt'))
Z3 = transpose(loadtxt('out/sol3.txt'))
Z4 = transpose(loadtxt('out/sol4.txt'))
Z5 = transpose(loadtxt('out/sol5.txt'))
Z6 = transpose(loadtxt('out/sol6.txt'))
Z7 = transpose(loadtxt('out/sol7.txt'))
abc = loadtxt('out/otherinfo.txt')
chg = loadtxt('out/chargeXYZ.txt')
col = loadtxt('out/chargeMag.txt')
print(abc)
print(chg)
ells = [Ellipse(xy=[0,0], width=2*abc[0], height=2*abc[1], angle=0) for i in range(8)]
area = 50

fig = plt.figure(0)

ellAlpha = .3
levs = [-5, -4, -2, -1, -.75, -.5, -.25, -.15, -.1, -.075, -.05, -.025, 0, .025, .05, .075, .1, .15, .2, .25, .5, .75, 1]
color = rnd.rand(3)
print('Color:')
print(color)

ax = fig.add_subplot(421, aspect='equal')
ax.add_artist(ells[0])
ells[0].set_clip_box(ax.bbox)
ells[0].set_alpha(ellAlpha)
ells[0].set_facecolor(color)
ax.contour(X, Y, Z0, levels=levs)
ax.scatter(chg[:,0],chg[:,1], s=area, c=col, alpha=0.5)
#ax.scatter(chg[0],chg[1], s=area, c=col, alpha=0.5)
ax.set_title('plot 0')


ax = fig.add_subplot(422, aspect='equal')
ax.contour(X, Y, Z1, levels=levs)
ax.scatter(chg[:,0],chg[:,1], s=area, c=col, alpha=0.5)
#ax.scatter(chg[0],chg[1], s=area, c=col, alpha=0.5)
ax.add_artist(ells[1])
ells[1].set_clip_box(ax.bbox)
ells[1].set_alpha(ellAlpha)
ells[1].set_facecolor(color)
ax.set_title('plot 1')


ax = fig.add_subplot(423, aspect='equal')
ax.contour(X, Y, Z2, levels=levs)
ax.scatter(chg[:,0],chg[:,1], s=area, c=col, alpha=0.5)
#ax.scatter(chg[0],chg[1], s=area, c=col, alpha=0.5)
ax.add_artist(ells[2])
ells[2].set_clip_box(ax.bbox)
ells[2].set_alpha(ellAlpha)
ells[2].set_facecolor(color)
ax.set_title('plot 2')


ax = fig.add_subplot(424, aspect='equal')
ax.contour(X, Y, Z3, levels=levs)
ax.scatter(chg[:,0],chg[:,1], s=area, c=col, alpha=0.5)
#ax.scatter(chg[0],chg[1], s=area, c=col, alpha=0.5)
ax.add_artist(ells[3])
ells[3].set_alpha(ellAlpha)
ells[3].set_facecolor(color)
ax.set_title('plot 3')


ax = fig.add_subplot(425, aspect='equal')
ax.contour(X, Y, Z4, levels=levs)
ax.scatter(chg[:,0],chg[:,1], s=area, c=col, alpha=0.5)
#scatter(chg[0],chg[1], s=area, c=col, alpha=0.5)
ax.add_artist(ells[4])
ells[4].set_alpha(ellAlpha)
ells[4].set_facecolor(color)
ax.set_title('plot 4')


ax = fig.add_subplot(426, aspect='equal')
ax.contour(X, Y, Z5, levels=levs)
ax.scatter(chg[:,0],chg[:,1], s=area, c=col, alpha=0.5)
#scatter(chg[0],chg[1], s=area, c=col, alpha=0.5)
ax.add_artist(ells[5])
ells[5].set_alpha(ellAlpha)
ells[5].set_facecolor(color)
ax.set_title('plot 5')


ax = fig.add_subplot(427, aspect='equal')
ax.contour(X, Y, Z6, levels=levs)
ax.scatter(chg[:,0],chg[:,1], s=area, c=col, alpha=0.5)
#ax.scatter(chg[0],chg[1], s=area, c=col, alpha=0.5)
ax.add_artist(ells[6])
ells[6].set_alpha(ellAlpha)
ells[6].set_facecolor(color)
ax.set_title('plot 6')

ax = fig.add_subplot(428, aspect='equal')
ax.contour(X, Y, Z7, levels=levs)
ax.scatter(chg[:,0],chg[:,1], s=area, c=col, alpha=0.5)
#ax.scatter(chg[0],chg[1], s=area, c=col, alpha=0.5)
ax.add_artist(ells[7])
ells[6].set_alpha(ellAlpha)
ells[6].set_facecolor(color)
ax.set_title('plot 7')

plt.show()
