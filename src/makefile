CLINKER=gcc

CFLAGS = -I. -lm -lgmp -lmpfr -llapack -lblas -lgsl -lgslcblas

SPH_OBJS = sphere/sphere.o sphere/sphSolv.o

ELL_OBJS = ellipsoid/ellipsoid.o ellipsoid/ellSolv.o


ex1: ${SPH_OBJS} ${ELL_OBJS} testProblem.o
	${CLINKER} -o $@.out $^ ${PETSC_LIB}
	${DSYMUTIL} $@

testell: ${ELL_OBJS} testEllipsoid.o
	${CLINKER} -o $@.out $^ ${PETSC_LIB}
	${DSYMUTIL} $@

testdassios: ${ELL_OBJS} testDassios.o
	${CLINKER} -o $@.out $^ ${PETSC_LIB}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
