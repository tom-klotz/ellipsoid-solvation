from numpy import *
from matplotlib.pyplot import *

figure()
A = loadtxt('plotsphtest4.txt')
for i in range(0,A.size):
    if (A.flat[i] > 1):
        A.flat[i] = 1;
print(A.max())
C = contour(A)
print(A.min())
colorbar(C, shrink=0.8, extend='both')
savefig('sphcplot.eps')
show()
