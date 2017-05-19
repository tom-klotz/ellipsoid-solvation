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

  mpfr_init(temp1);
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


#undef __FUNCT__
#define __FUNCT__ "VecViewMPFR"
PetscErrorCode VecViewMPFR(PetscInt size_d, mpfr_t *vec) {
  PetscErrorCode ierr;
  double val;
  PetscFunctionBegin;
  printf("-------------------\n");
  for(PetscInt k=0; k < size_d; ++k) {
    val = mpfr_get_d(vec[k], MPFR_RNDN);
    printf("vec[%d] = %4.4f\n", k, val);
  }
  printf("-------------------\n");
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatViewMPFR"
PetscErrorCode MatViewMPFR(PetscInt m, PetscInt n, mpfr_t *mat) {
  PetscErrorCode ierr;
  double val;
  PetscFunctionBegin;
  printf("-------------------\n");
  for(PetscInt i=0; i < m; ++i) {
    for(PetscInt j=0; j < n; ++j) {
      val = mpfr_get_d(mat[i*n+j], MPFR_RNDN);
      printf("%3.3f ", val);
    }
    printf("\n");
  }
  printf("-------------------\n");

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Vec2NormMPFR"
PetscErrorCode Vec2NormMPFR(PetscInt size_d, mpfr_t *vec, mpfr_t *norm)
{
  PetscErrorCode ierr;
  mpfr_t val;
  PetscFunctionBegin;
  mpfr_init(val);
  mpfr_set_d(*norm, 0.0, MPFR_RNDN);
  for(PetscInt k=0; k < size_d; ++k) {
    mpfr_mul(val, vec[k], vec[k], MPFR_RNDN);
    mpfr_add(*norm, *norm, val, MPFR_RNDN);
  }
  mpfr_sqrt(*norm, *norm, MPFR_RNDN);

  
  mpfr_clear(val);
  PetscFunctionReturn(0);
}


//returns (lambda/|v1|^2)*v1*v1^T*vec
#undef __FUNCT__
#define __FUNCT__ "OrthogMap"
PetscErrorCode OrthogMap(PetscInt size_d, mpfr_t *x, mpfr_t *v1, mpfr_t *lambda,  mpfr_t *result) {

  PetscErrorCode ierr;
  mpfr_t temp1, temp2, coef;
  PetscFunctionBegin;

  mpfr_inits(temp1, temp2, coef, NULL);
  
  ierr = Vec2NormMPFR(size_d, v1, &temp1);CHKERRQ(ierr);
  mpfr_mul(temp1, temp1, temp1, MPFR_RNDN);
  mpfr_div(coef, *lambda, temp1, MPFR_RNDN);
  
  
  //compute v1^T*vec
  mpfr_set_d(temp2, 0.0, MPFR_RNDN);
  for(PetscInt k=0; k < size_d; ++k) {
    mpfr_mul(temp1, v1[k], x[k], MPFR_RNDN);
    mpfr_add(temp2, temp2, temp1, MPFR_RNDN);
  }
  mpfr_mul(coef, coef, temp2, MPFR_RNDN);

  //multiply v1 by coef and store to result
  for(PetscInt k=0; k < size_d; ++k)
    mpfr_mul(result[k], v1[k], coef, MPFR_RNDN);
  
  mpfr_clears(temp1, temp2, coef, NULL);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "eigsMPFR"
PetscErrorCode eigsMPFR(PetscInt size_d, mpfr_t *mat, mpfr_t *result) {

  PetscErrorCode ierr;
  mpfr_t temp1;
  mpfr_t norm;
  mpfr_t error, error2;
  mpfr_t *tempvec, *tempvec2;
  mpfr_t *eigvec;
  mpfr_t *prev;
  mpfr_t *eigvecs;
  mpfr_t *eigvals;
  PetscFunctionBegin;

  mpfr_inits(norm, error, error2, NULL);
  
  eigvec = (mpfr_t*) malloc(sizeof(mpfr_t)*size_d);
  prev = (mpfr_t*) malloc(sizeof(mpfr_t)*size_d);
  eigvals = (mpfr_t*) malloc(sizeof(mpfr_t)*size_d);
  tempvec = (mpfr_t*) malloc(sizeof(mpfr_t)*size_d);
  tempvec2 = (mpfr_t*) malloc(sizeof(mpfr_t)*size_d);
  
  
  for(PetscInt k=0; k < size_d; ++k) {
    mpfr_inits(tempvec2[k], tempvec[k], eigvals[k], eigvec[k], prev[k], NULL);
  }
  for(PetscInt i=0; i < size_d; ++i) {
    for(PetscInt l=0; l < size_d; ++l) {
      mpfr_set_d(eigvec[l], 0.5, MPFR_RNDN);
    }
    mpfr_set_d(error, 1.0, MPFR_RNDN);
    mpfr_set_d(error2, 1.0, MPFR_RNDN);
    mpfr_log10(error, error, MPFR_RNDN);
    mpfr_log10(error2, error2, MPFR_RNDN);
    PetscInt j=0;
    while(mpfr_get_d(error, MPFR_RNDN) > -32 && mpfr_get_d(error2, MPFR_RNDN) > -32) {
      j++;
      for(PetscInt l=0; l < size_d; ++l)
	mpfr_set(prev[l], eigvec[l], MPFR_RNDN);
      ierr = MatVecMPFR(size_d, mat, eigvec, tempvec2);CHKERRQ(ierr);
      //subtract off already found eigenvectors
      for(PetscInt k=0; k<i; ++k) {
	//PetscErrorCode OrthogMap(PetscInt size_d, mpfr_t *x, mpfr_t *v1, mpfr_t *lambda,  mpfr_t *result) {
	OrthogMap(size_d, eigvec, result+k*size_d, eigvals+k, tempvec);
	for(PetscInt l=0; l<size_d; ++l) {
	  mpfr_sub(tempvec2[l], tempvec2[l], tempvec[l], MPFR_RNDN);
	}
      }
      //normalize eig vector
      ierr = Vec2NormMPFR(size_d, tempvec2, &norm);CHKERRQ(ierr);
      for(PetscInt k=0; k < size_d; ++k)
	mpfr_div(eigvec[k], tempvec2[k], norm, MPFR_RNDN);
      //compute error
      for(PetscInt k=0; k < size_d; ++k) {
	mpfr_sub(tempvec[k], eigvec[k], prev[k], MPFR_RNDN);
	mpfr_add(tempvec2[k], eigvec[k], prev[k], MPFR_RNDN);
      }
      ierr = Vec2NormMPFR(size_d, tempvec, &error);CHKERRQ(ierr);
      ierr = Vec2NormMPFR(size_d, tempvec2, &error2);CHKERRQ(ierr);
      ierr = mpfr_log10(error, error, MPFR_RNDN);
      ierr = mpfr_log10(error2, error2, MPFR_RNDN);
      //if(i==1) {
      //printf("e: %2.2f\n", mpfr_get_d(error, MPFR_RNDN));
      //printf("norm: %3.3e\n", mpfr_get_d(norm, MPFR_RNDN));
      //ierr = VecViewMPFR(size_d, eigvec);CHKERRQ(ierr);
      //}
     
	
    }
    //if(i==0) {
    //ierr = VecViewMPFR(size_d, eigvec);CHKERRQ(ierr);
    //}
    //add newly found eigvec to eigvecss
    for(PetscInt k=0; k < size_d; ++k)
      mpfr_set(result[i*size_d+k], eigvec[k], MPFR_RNDN);
    //add new eigval
    mpfr_set(eigvals[i], norm, MPFR_RNDN);
    
    printf("EIGVAL IS %3.3f\n", mpfr_get_d(norm, MPFR_RNDN));

  }
  mpfr_clear(norm);
  PetscFunctionReturn(0);
}

