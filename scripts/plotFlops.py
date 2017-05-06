from numpy import *
from matplotlib.pyplot import *

figure()
A = loadtxt('flopsext.txt', usecols={3})
B = loadtxt('flopsint.txt', usecols={3})
A = cumsum(A)
Awow = A / 100;
B = cumsum(B)
semilogy(A, '+', Awow, '_', B, 'o')
show()
