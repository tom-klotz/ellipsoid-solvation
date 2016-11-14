#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpfr.h>
#include <petsc.h>

#include "ellipsoid/ellipsoid.h"
#include "tanhsinh.h"


#undef __FUNCT__
#define __FUNCT__ "NormConstantIntFixed"
PetscErrorCode NormConstantIntFixed(EllipsoidalSystem *e, PetscInt n, PetscInt p, PetscInt prec, PetscInt nPts, PetscReal *normConst)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;


  FuncInfo2 ctx1 = { .e = e, .n = n, .p = p, .numeratorType = 0, .denomSign = 1 };
  FuncInfo2 ctx2 = { .e = e, .n = n, .p = p, .numeratorType = 1, .denomSign = 1 };
  FuncInfo2 ctx3 = { .e = e, .n = n, .p = p, .numeratorType = 0, .denomSign = -1 };
  FuncInfo2 ctx4 = { .e = e, .n = n, .p = p, .numeratorType = 1, .denomSign = -1 };

  PetscReal integrals[4];

  mpfr_t mpfrzero;
  mpfr_t mpfrone;
  mpfr_inits(mpfrzero, mpfrone, NULL);
  mpfr_set_d(mpfrzero, 0.0, MPFR_RNDN);
  mpfr_set_d(mpfrone, 1.0, MPFR_RNDN);
  
  ierr = DEQuad((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, e->hp_h, e->hp_k, prec, nPts, integrals+0, &ctx1);

  ierr = DEQuad((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, e->hp_h, e->hp_k, prec, nPts, integrals+1, &ctx2);

  ierr = DEQuad((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, mpfrzero, e->hp_h, prec, nPts, integrals+2, &ctx3);

  ierr = DEQuad((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, mpfrzero, e->hp_h, prec, nPts, integrals+3, &ctx4);

  *normConst = 8.0*(integrals[2]*integrals[1] - integrals[0]*integrals[3]);

  mpfr_clears(mpfrzero, mpfrone, NULL);
  
  ierr = PetscLogFlops(4);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "NormPlot"
PetscErrorCode NormPlot()
{
  const PetscInt NUM_SOLUTIONS = 20;
  const PetscInt POINTS_MIN = 5;
  const PetscInt POINTS_STEP = 5;
  const PetscInt prec = 16;
  const PetscReal a = 3.0;
  const PetscReal b = 2.0;
  const PetscReal c = 1.0;
  PetscReal sol1, sol2, sol3, sol4, solOld;
  EllipsoidalSystem e;
  PetscErrorCode ierr;

  initEllipsoidalSystem(&e, a, b, c);

  ierr = NormConstantIntFixed(&e, 3, 2, prec, 10, &sol1);CHKERRQ(ierr);
  printf("new norm constant: %15.15f\n", sol1);
  ierr = NormConstantIntFixed(&e, 3, 2, prec, 20, &sol2);CHKERRQ(ierr);
  printf("new norm constant: %15.15f\n", sol2);
  ierr = NormConstantIntFixed(&e, 3, 2, prec, 40, &sol3);CHKERRQ(ierr);
  printf("new norm constant: %15.15f\n", sol3);
  ierr = NormConstantIntFixed(&e, 3, 2, prec, 80, &sol4);CHKERRQ(ierr);
  printf("new norm constant: %15.15f\n", sol4);
  ierr = calcNormalization(&e, 3, 2, &solOld);CHKERRQ(ierr);
  printf("old norm constant: %15.15f\n", solOld);
  
  PetscFunctionBegin;
  



  
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "main"
PetscErrorCode main(int argc, char **argv)
{

  const PetscInt NUM_DIGITS = 64;
  const PetscInt PRECISION  = 16;
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

  /*
  int fac = 1;
  for(int i=0; i<10; ++i) {
    printf("\n########################################################\n");
    ierr = DEQuad(inte, a, b, PRECISION, fac*nPts, &val, NULL);CHKERRQ(ierr);
    fac *= 2;
  }
  */
  ierr = NormPlot();CHKERRQ(ierr);
  
  ierr = PetscFinalize();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
