#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include <petsc.h>
#include "ellipsoid.h"
#include "../sphere/sphere.h"
#include "ellSolv.h"


#undef __FUNCT__
#define __FUNCT__ "calcSolidHarmonic"
PetscErrorCode calcSolidHarmonic(EllipsoidalSystem* e, PetscReal lambda, PetscReal mu, PetscReal nu, PetscInt n, PetscInt p, PetscReal *val)
{
  PetscInt signm = 1;
  PetscInt signn = 1;
  PetscReal eL, eM, eN;
  PetscFunctionBegin;
  
  if(mu < 0)
    signm = -1;
  if(nu < 0)
    signn = -1;

  eL = calcLame(e, n, p, lambda, signm, signn);
  eM = calcLame(e, n, p, mu, signm, signn);
  eN = calcLame(e, n, p, nu, signm, signn);
  
  *val = eL*eM*eN;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "calcSolidHarmonicVec"
/*
  ellPoints is ellPoints of size 3*nPoints containing ellipsoidal coordinates
  outputs solid harmonics to values of size nPoints
*/
PetscErrorCode calcSolidHarmonicVec(EllipsoidalSystem* e, PetscInt nPoints, Vec ellPoints, PetscInt n, PetscInt p, Vec values)
{
  PetscErrorCode ierr;
  PetscInt signm = 1;
  PetscInt signn = 1;
  PetscReal lambda, mu, nu;
  PetscReal eL, eM, eN;
  PetscReal Enp;
  const PetscScalar* ellPointsArray;
  PetscFunctionBegin;

  ierr = VecGetArrayRead(ellPoints, &ellPointsArray); CHKERRQ(ierr);
  for(PetscInt k=0; k<nPoints; ++k) {
    lambda = ellPointsArray[3*k+0];
    mu     = ellPointsArray[3*k+1];
    nu     = ellPointsArray[3*k+2];
    if(mu < 0)
      signm = -1;
    else signm = 1;
    if(nu < 0)
      signn = -1;
    else signn = 1;
    eL = calcLame(e, n, p, lambda, signm, signn);
    eM = calcLame(e, n, p, mu, signm, signn);
    eN = calcLame(e, n, p, nu, signm, signn);
    Enp = eL*eM*eN;
    ierr = VecSetValue(values, k, Enp, INSERT_VALUES); CHKERRQ(ierr);
  }
  ierr = VecAssemblyBegin(values); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(values); CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(ellPoints, &ellPointsArray); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "calcCoulombEllCoefs"
PetscErrorCode calcCoulombEllCoefs(EllipsoidalSystem* e, PetscInt nCharges, Vec srcPoints, Vec srcCharges, PetscInt Nmax, Vec coefs)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  
  PetscFunctionReturn(0);
}
