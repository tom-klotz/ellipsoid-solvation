CLINKER=gcc

CFLAGS   = -I../include -I. -std=c99
SPH_LIBS = -lm -lgmp -lmpfr -llapack -lblas -L../lib -lgsl -lgslcblas

SPH_OBJS = ,,/sphere/sphere.o sphere/sphSolv.o

ELL_OBJS = ../ellipsoid/ellipsoid.o ../ellipsoid/ellSolv.o ../ellipsoid/ellSolvPetsc.o ../ellipsoid/MPFReigs.o


testEllMPFR: ${ELL_OBJS} testEllMPFR.o
	${CLINKER} -o $@ $^ ${SPH_LIBS} ${PETSC_LIB}
	${DSYMUTIL} $@


include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
