all:
	gcc -O3 -g -o main ellipsoid.c sphere.c ellSolv.c sphSolv.c testProblem.c -lm -lgmp -lmpfr -llapack -lblas -lgsl -lgslcblas
test:
	gcc -O3 -g -o main ellipsoid.c testProblem.c -lm -lgmp -lmpfr -llapack -lblas -lgsl -lgslcblas
speed:
	gcc -O3 -g -o main ellipsoid.c speedtester.c -lgmp -lmpfr -llapack -lblas
testell:
	gcc -O3 -g -o main ellipsoid.c testEllipsoid.c -lm -lgmp -lmpfr -llapack -lblas
testdassios:
	gcc -O3 -g -o main ellipsoid.c testDassios.c -lm -lgmp -lmpfr -llapack -lblas
clean:
	rm main
