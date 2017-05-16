#include <math.h>
#include <stdio.h>
#include <mpfr.h>
#include <stdlib.h>
#include <string.h>
#include <petsc.h>
#include <time.h>


#undef __FUNCT__
#define __FUNCT__ "MatVecMPFR"
PetscErrorCode MatVecMPFR(PetscInt size_d, mpfr_t *Mat, mpfr_t *vec, mpfr_t *result) {

  PetscErrorCode ierr;
  mpfr_t temp1;
  mpfr_t *res;
  PetscFunctionBegin;

  res = (mpfr_t*) malloc(sizeof(mpfr_t)*size_d);
  
  for(PetscInt k=0; k < size_d; ++k)
    mpfr_init(res[k]);
  
  for(PetscInt i=0; i<size_d; ++i) {
    mpfr_set_d(res[i], 0.0, MPFR_RNDN);
    for(PetscInt j=0; j<size_d; ++j) {
      mpfr_mul(temp1, Mat[i*size_d+j], vec[j], MPFR_RNDN);
      mpfr_add(res[i], res[i], temp1, MPFR_RNDN);
    }
  }
  for(PetscInt k=0; k < size_d; ++k) {
    mpfr_set(result[k], res[k], MPFR_RNDN);
    mpfr_clear(res[k]);
  }
  
  mpfr_clear(temp1);
  free(res);
  PetscFunctionReturn(0);
}

