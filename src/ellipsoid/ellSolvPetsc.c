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
#define __FUNCT__ "CalcSolidHarmonic"
PetscErrorCode CalcSolidHarmonic(EllipsoidalSystem* e, PetscReal lambda, PetscReal mu, PetscReal nu, PetscInt n, PetscInt p, PetscReal *val)
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
#define __FUNCT__ "CalcSolidHarmonicVec"
/*
  ellPoints is ellPoints of size 3*nPoints containing ellipsoidal coordinates
  outputs solid harmonics to values of size nPoints
*/
PetscErrorCode CalcSolidHarmonicVec(EllipsoidalSystem* e, PetscInt nPoints, Vec ellPoints, PetscInt n, PetscInt p, Vec values)
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
#define __FUNCT__ "CalcCoulombEllCoefs"
/*
  GnpVals should be uninitialized. I
*/
PetscErrorCode CalcCoulombEllCoefs(EllipsoidalSystem* e, PetscInt nPoints, Vec srcPoints, Vec srcCharges, PetscInt Nmax, Vec *GnpVals)
{
  PetscErrorCode ierr;
  PetscReal normConstant;
  PetscReal Gnp;
  PetscInt count;

  Vec EnpVals;
  const PetscScalar* EnpValsArray;
  const PetscScalar* srcPointsArray;
  const PetscScalar* srcChargesArray;
  PetscFunctionBegin;

  
  ierr = VecCreateSeq(PETSC_COMM_SELF, nPoints, &EnpVals);CHKERRQ(ierr);

  count = 0;
  for(PetscInt n=0; n<=Nmax; ++n) {
    for(PetscInt p=0; p<2*n+1; ++p)
      count++;
  }
  ierr = VecDestroy(GnpVals);CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, count, GnpVals);CHKERRQ(ierr);
  
  ierr = VecGetArrayRead(srcCharges, &srcChargesArray);CHKERRQ(ierr);
  count = 0;
  for(PetscInt n=0; n<=Nmax; ++n) {
    for(PetscInt p=0; p<2*n+1; ++p) {
      normConstant = calcNormalization(e, n, p);
      
      ierr = CalcSolidHarmonicVec(e, nPoints, srcPoints, n, p, EnpVals);
      
      ierr = VecGetArrayRead(srcPoints, &srcPointsArray);CHKERRQ(ierr);
      ierr = VecGetArrayRead(EnpVals, &EnpValsArray);CHKERRQ(ierr);
      Gnp = 0;
      for(PetscInt k=0; k<nPoints; ++k) {
	Gnp += srcChargesArray[k]*EnpValsArray[k];
      }
      ierr = VecRestoreArrayRead(EnpVals, &EnpValsArray);CHKERRQ(ierr);
      ierr = VecRestoreArrayRead(srcPoints, &srcPointsArray);CHKERRQ(ierr);
      Gnp *= (4*M_PI)/((2.0*n+1.0)*normConstant);
      ierr = VecSetValues(*GnpVals, 1, &count, &Gnp, INSERT_VALUES);CHKERRQ(ierr);
      count++;
    }
  }
  ierr = VecAssemblyBegin(*GnpVals);CHKERRQ(ierr); ierr = VecAssemblyEnd(*GnpVals);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(srcCharges, &srcChargesArray);CHKERRQ(ierr);

  ierr = VecDestroy(&EnpVals); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
