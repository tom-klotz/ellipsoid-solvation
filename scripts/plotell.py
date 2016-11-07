from numpy import *
from matplotlib.pyplot import *

figure()
A = loadtxt('solvPlot.txt', skiprows=2)
#for i in range(0,A.size):
    #if (A.flat[i] > 1):
    #A.flat[i] = 1;
print(A)
C = contour(A)
print(A.min())
colorbar(C, shrink=0.8, extend='both')
savefig('ellcplot.eps')
show()
