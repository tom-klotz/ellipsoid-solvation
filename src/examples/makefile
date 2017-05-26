CLINKER=gcc

CFLAGS   = -I../include -I. -std=c99
SPH_LIBS = -lm -lgmp -lmpfr -llapack -lblas -L../lib -lgsl -lgslcblas

QUAD_OBJS = ../tanhsinh.o

SPH_OBJS = ../sphere/sphere.o ../sphere/sphSolv.o

ELL_OBJS = ../ellipsoid/ellipsoid.o ../ellipsoid/ellSolv.o ../ellipsoid/ellSolvPetsc.o ../ellipsoid/MPFReigs.o


ex1: ${SPH_OBJS} ${ELL_OBJS} ex1.o testFunctions.o
	${CLINKER} -o $@ $^ ${SPH_LIBS} ${PETSC_LIB}
	${DSYMUTIL} $@

ex2: ${ELL_OBJS} ex2.o
	${CLINKER} -o $@ $^ ${SPH_LIBS} ${PETSC_LIB}
	${DSYMUTIL} $@

ex3: ${ELL_OBJS} ${QUAD_OBJS} ex3.o
	${CLINKER} -o $@ $^ ${SPH_LIBS} ${PETSC_LIB}
	${DSYMUTIL} $@


include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
