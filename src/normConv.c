#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpfr.h>
#include <petsc.h>

#include "ellipsoid/ellipsoid.h"
#include "tanhsinh.h"


#undef __FUNCT__
#define __FUNCT__ "main"
PetscErrorCode main(int argc, char **argv)
{

  const PetscInt NUM_DIGITS = 64;
  const PetscInt PRECISION  = 32;
  mpfr_t a, b;
  mpfr_t val;
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  const PetscInt nPts = 10;

  /* initializations */
  ierr = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRQ(ierr);
  ierr = PetscLogDefaultBegin();CHKERRQ(ierr);
  
  mpfr_set_default_prec(4*NUM_DIGITS);
  mpfr_inits(a, b, val, NULL);
  mpfr_set_d(a, 2.0, MPFR_RNDN);
  mpfr_set_d(b, 3.0, MPFR_RNDN);


  
  ierr = DEQuad(inte, a, b, PRECISION, nPts, &val, NULL);CHKERRQ(ierr);
  printf("\n########################################################\n");
  ierr = DEQuad(inte, a, b, PRECISION, 2*nPts, &val, NULL);CHKERRQ(ierr);
  
  ierr = PetscFinalize();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
