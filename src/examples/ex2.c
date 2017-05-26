#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpfr.h>
#include <petsc.h>

#include "../ellipsoid/ellipsoid.h"
#include "../ellipsoid/ellSolv.h"
#include "../constants.h"
#include "testFunctions.h"



#undef __FUNCT__
#define __FUNCT__ "main"
PetscErrorCode main(int argc, char **argv)
{

  PetscErrorCode ierr;
  EllipsoidalSystem e;
  PetscInt prec = 64;
  PetscReal a = 3.0;
  PetscReal b = 2.0;
  PetscReal c = 1.0;
  double doubleEnp;
  mpfr_t x;
  mpfr_t mpfrEnp;
  PetscFunctionBeginUser;

  ierr = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRQ(ierr);
  ierr = PetscLogDefaultBegin();CHKERRQ(ierr);

  
  mpfr_inits(mpfrEnp, x, NULL);
  
  ierr = initEllipsoidalSystem(&e, a, b, c, prec);CHKERRQ(ierr);

  mpfr_set_d(x, 1.0, MPFR_RNDN);
  //calculate with double
  ierr = calcLame(&e, 3, 4, 1.0, 1, 1, &doubleEnp);
  //calculate with high precision
  ierr = CalcLameMPFR(&e, 3, 4, x, 1, 1, &mpfrEnp);

  printf("Enp(double): %18.18f\n", doubleEnp);
  printf("Enp(mpfr):   %18.18f\n", mpfr_get_d(mpfrEnp, MPFR_RNDN));

  mpfr_clears(mpfrEnp, x, NULL);
  ierr = PetscFinalize();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
