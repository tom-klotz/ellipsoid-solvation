all:
	g++ -O3 -g -o main ellipsoid.c testProblem.cpp -lgmp -lmpfr -llapack -lblas
test:
	g++ -O3 -g -o main ellipsoid.c testProblem.c -lgmp -lmpfr -llapack -lblas
testcpp:
	g++ -O3 -g -o main ellipsoid.c testProblem.cpp -lgmp -lmpfr -llapack -lblas
speed:
	g++ -O3 -g -o main ellipsoid.c speedtester.c -lgmp -lmpfr -llapack -lblas
ex1:
	g++ -o main ellipsoid.c exteriorSolve.cpp -lgmp -lmpfr -llapack -lblas -lboost_system
testell:
	g++ -O2 -g -o main ellipsoid.c testEllipsoid.c -lgmp -lmpfr -llapack -lblas
testdassios:
	g++ -g -o main ellipsoid.c testDassios.c -lgmp -lmpfr -llapack -lblas
clean:
	rm main