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
#define __FUNCT__ "CalcEllipsoidSolvationPotential"
PetscErrorCode CalcEllipsoidSolvationPotential(PetscReal a, PetscReal b, PetscReal c, PetscReal eps1, PetscReal eps2, PetscInt nSource, Vec sourceXYZ, Vec sourceMag, PetscInt nTarget, Vec targetXYZ, PetscInt Nmax, Vec targetSol)
{
  PetscErrorCode ierr;
  EllipsoidalSystem e;
  PetscReal x, y, z;
  PetscInt n, p;
  PetscInt intPts, extPts;
  Vec tarIntXYZ, tarExtXYZ;
  PetscScalar *tarIntXYZArray, *tarExtXYZArray;
  Vec tarIntEll, tarExtEll;
  Vec srcEll, tarEll;
  Vec coulCoefs, reactCoefs, extCoefs;
  Vec EnpVals, FnpVals;
  PetscInt *isExt;
  const PetscScalar *targetXYZArray, *coulCoefsArray, *reactCoefsArray, *extCoefsArray;
  PetscScalar *targetSolArray;
  const PetscScalar *EnpValsArray, *FnpValsArray;

  PetscLogStage stageSolids;
  PetscLogEvent INTERIOR_FLOPS;
  PetscLogEvent EXTERIOR_FLOPS;
  PetscLogDouble user_log_flops;
  PetscFunctionBegin;

  ierr = PetscLogStageRegister("Solid Harmonic Flops", &stageSolids);
  
  initEllipsoidalSystem(&e, a, b, c);

  // calculate the number of interior and exterior points
  ierr = VecGetArrayRead(targetXYZ, &targetXYZArray);CHKERRQ(ierr);
  intPts = 0; extPts = 0;
  for(int k=0; k < nTarget; ++k) {
    x = targetXYZArray[3*k+0];
    y = targetXYZArray[3*k+1];
    z = targetXYZArray[3*k+2];
    if((x*x)/(a*a) + (y*y)/(b*b) + (z*z)/(c*c) <= 1)
      intPts++;
    else
      extPts++;
  }
  // init interior and exterior target point vectors
  ierr = VecCreateSeq(PETSC_COMM_SELF, 3*intPts, &tarIntXYZ);CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, 3*extPts, &tarExtXYZ);CHKERRQ(ierr);
  // sort points into int and ext vectors
  ierr = VecGetArray(tarIntXYZ, &tarIntXYZArray);CHKERRQ(ierr);
  ierr = VecGetArray(tarExtXYZ, &tarExtXYZArray);CHKERRQ(ierr);

  ierr = PetscMalloc1(sizeof(PetscInt)*nTarget, &isExt);CHKERRQ(ierr);
  
  intPts = 0; extPts = 0;
  for(int k=0; k < nTarget; ++k) {
    x = targetXYZArray[3*k+0];
    y = targetXYZArray[3*k+1];
    z = targetXYZArray[3*k+2];
    if((x*x)/(a*a) + (y*y)/(b*b) + (z*z)/(c*c) <= 1) {
      tarIntXYZArray[3*intPts+0] = x;
      tarIntXYZArray[3*intPts+1] = y;
      tarIntXYZArray[3*intPts+2] = z;
      isExt[k] = 0;
      intPts++;   
    }
    else {
      tarExtXYZArray[3*extPts+0] = x;
      tarExtXYZArray[3*extPts+1] = y;
      tarExtXYZArray[3*extPts+2] = z;
      isExt[k] = 1;
      extPts++;
    }
  }
  ierr = VecRestoreArray(tarIntXYZ, &tarIntXYZArray);CHKERRQ(ierr);
  ierr = VecRestoreArray(tarExtXYZ, &tarExtXYZArray);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(targetXYZ, &targetXYZArray);CHKERRQ(ierr);
  
  // create source ellipsoidal vec and convert from xyz
  ierr = VecCreateSeq(PETSC_COMM_SELF, 3*nSource, &srcEll);CHKERRQ(ierr);
  ierr = CartesianToEllipsoidalVec(&e, sourceXYZ, srcEll);CHKERRQ(ierr);
  // create target ellipsoidal vec and convert from xyz
  ierr = VecCreateSeq(PETSC_COMM_SELF, 3*nTarget, &tarEll);CHKERRQ(ierr);
  ierr = CartesianToEllipsoidalVec(&e, targetXYZ, tarEll);CHKERRQ(ierr);
  // create target interior ellipsoidal vec and convert from xyz
  ierr = VecCreateSeq(PETSC_COMM_SELF, 3*intPts, &tarIntEll);CHKERRQ(ierr);
  ierr = CartesianToEllipsoidalVec(&e, tarIntXYZ, tarIntEll);CHKERRQ(ierr);
  // create target exterior ellipsoidal vec and convert from xyz
  ierr = VecCreateSeq(PETSC_COMM_SELF, 3*extPts, &tarExtEll);CHKERRQ(ierr);
  ierr = CartesianToEllipsoidalVec(&e, tarExtXYZ, tarExtEll);CHKERRQ(ierr);
  
  // calculate coulomb coefficients
  ierr = CalcCoulombEllCoefs(&e, nSource, srcEll, sourceMag, Nmax, &coulCoefs);CHKERRQ(ierr);
  // init reacCoefs, extCoefs from size of coulCoefs
  ierr = VecDuplicate(coulCoefs, &reactCoefs);CHKERRQ(ierr);
  ierr = VecDuplicate(coulCoefs, &extCoefs);CHKERRQ(ierr);
  // calc reaction and exterior expansion coefs from coulomb coefs
  ierr = CalcReactAndExtCoefsFromCoulomb(&e, eps1, eps2, Nmax, coulCoefs, reactCoefs, extCoefs);


  // create EnpVals and FnpVals vectors
  ierr = VecCreateSeq(PETSC_COMM_SELF, intPts, &EnpVals);CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, extPts, &FnpVals);CHKERRQ(ierr);

  
  // get read-only pointers for expansion coefficients
  ierr = VecGetArrayRead(coulCoefs, &coulCoefsArray);CHKERRQ(ierr);
  ierr = VecGetArrayRead(reactCoefs, &reactCoefsArray);CHKERRQ(ierr);
  ierr = VecGetArrayRead(extCoefs, &extCoefsArray);CHKERRQ(ierr);
  // get write pointer for solution vector
  ierr = VecZeroEntries(targetSol);CHKERRQ(ierr);
  ierr = VecGetArray(targetSol, &targetSolArray);CHKERRQ(ierr);
  // loop
  PetscInt ind = 0;
  PetscInt intInd;
  PetscInt extInd;
  ierr = PetscLogEventRegister("Interior Harmonic Flops", 0, &INTERIOR_FLOPS);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("Exterior Harmonic Flops", 0, &EXTERIOR_FLOPS);CHKERRQ(ierr);
  ierr = PetscLogEventBegin(INTERIOR_FLOPS, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = PetscLogEventBegin(EXTERIOR_FLOPS, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = PetscLogEventDeactivate(INTERIOR_FLOPS);CHKERRQ(ierr);
  ierr = PetscLogEventDeactivate(EXTERIOR_FLOPS);CHKERRQ(ierr);
  for(n=0; n <= Nmax; ++n) {
    printf("%d/%d\n", n, Nmax);
    for(p=0; p < 2*n+1; ++p) {
      
      PetscReal Gnp =  coulCoefsArray[ind];
      PetscReal Cnp =   extCoefsArray[ind];
      PetscReal Bnp = reactCoefsArray[ind];
      // calc solid harmonics for interior,exterior points


      ierr = PetscLogStagePush(stageSolids);CHKERRQ(ierr);
      ierr = PetscLogEventBegin(INTERIOR_FLOPS,0,0,0,0);CHKERRQ(ierr);
      ierr = CalcSolidInteriorHarmonicVec(&e, intPts, tarIntEll, n, p, EnpVals);CHKERRQ(ierr);
      ierr = PetscLogEventEnd  (INTERIOR_FLOPS,0,0,0,0);CHKERRQ(ierr);
      ierr = PetscLogEventBegin(EXTERIOR_FLOPS,0,0,0,0);CHKERRQ(ierr);
      ierr = CalcSolidExteriorHarmonicVec(&e, extPts, tarExtEll, n, p, FnpVals);CHKERRQ(ierr);
      ierr = PetscLogEventEnd(EXTERIOR_FLOPS,0,0,0,0);CHKERRQ(ierr);
      ierr = PetscLogStagePop();CHKERRQ(ierr);
      
      ierr = VecGetArrayRead(EnpVals, &EnpValsArray);CHKERRQ(ierr);
      ierr = VecGetArrayRead(FnpVals, &FnpValsArray);CHKERRQ(ierr);
      intInd = 0;
      extInd = 0;
      for(PetscInt k=0; k<nTarget; ++k) {

	if(isExt[k] == 0) { // interior points
	  //printf("E%d%d: %15.15f\nBnp %d%d: %15.15f\n\n", n, p, n, p, EnpValsArray[intInd], reactCoefsArray[ind]);
	  targetSolArray[k] += reactCoefsArray[ind]*EnpValsArray[intInd];
	  //ierr = CalcSolidInterior
	  //targetSolArray[k] += Gnp*EnpValsArray[intInd];
	  //targetSolArray[k] += 0;
	  intInd++;
	}
	else { // exterior points
	  //printf("F%d%d: %15.15f\nCnp %d%d: %15.15f\n\n", n, p, n, p, FnpValsArray[intInd], extCoefsArray[ind]);
	  targetSolArray[k] += extCoefsArray[ind]*FnpValsArray[extInd];
	  //targetSolArray[k] += (Gnp*FnpValsArray[extInd])/eps1;
	  targetSolArray[k] += 0;
	  extInd++;
	}
      }
      ierr = VecRestoreArrayRead(EnpVals, &EnpValsArray);CHKERRQ(ierr);
      ierr = VecRestoreArrayRead(FnpVals, &FnpValsArray);CHKERRQ(ierr);
      
      ind++;
      
    }
  }
  //ierr = PetscLogEventEnd(SOLID_FLOPS, 0, 0, 0, 0);

  ierr = VecRestoreArray(targetSol, &targetSolArray);CHKERRQ(ierr);
  // restore read-only pointers for expansion coefficients
  ierr = VecGetArrayRead(coulCoefs , &coulCoefsArray );CHKERRQ(ierr);
  ierr = VecGetArrayRead(reactCoefs, &reactCoefsArray);CHKERRQ(ierr);
  ierr = VecGetArrayRead(extCoefs  , &extCoefsArray  );CHKERRQ(ierr);




  ierr = VecDestroy(&tarIntXYZ);CHKERRQ(ierr);CHKERRQ(ierr);
  ierr = VecDestroy(&tarExtXYZ);CHKERRQ(ierr);CHKERRQ(ierr);
  ierr = VecDestroy(&tarIntEll);CHKERRQ(ierr);CHKERRQ(ierr);
  ierr = VecDestroy(&tarExtEll);CHKERRQ(ierr);CHKERRQ(ierr);
  ierr = VecDestroy(&srcEll);CHKERRQ(ierr);CHKERRQ(ierr);
  ierr = VecDestroy(&tarEll);CHKERRQ(ierr);CHKERRQ(ierr);
  ierr = VecDestroy(&coulCoefs);CHKERRQ(ierr);
  ierr = VecDestroy(&reactCoefs);CHKERRQ(ierr);
  ierr = VecDestroy(&extCoefs);CHKERRQ(ierr);
  ierr = VecDestroy(&EnpVals); CHKERRQ(ierr);
  ierr = VecDestroy(&FnpVals); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "CalcSolidInteriorHarmonic"
PetscErrorCode CalcSolidInteriorHarmonic(EllipsoidalSystem* e, PetscReal lambda, PetscReal mu, PetscReal nu, PetscInt n, PetscInt p, PetscReal *val)
{
  PetscErrorCode ierr;
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
  
  ierr= PetscLogFlops(3);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "CalcSolidInteriorHarmonicVec"
/*
  ellPoints is ellPoints of size 3*nPoints containing ellipsoidal coordinates
  outputs solid harmonics to values of size nPoints
*/
PetscErrorCode CalcSolidInteriorHarmonicVec(EllipsoidalSystem* e, PetscInt nPoints, Vec ellPoints, PetscInt n, PetscInt p, Vec values)
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
  
  ierr = PetscLogFlops(3*nPoints);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CalcSolidExteriorHarmonic"
PetscErrorCode CalcSolidExteriorHarmonic(EllipsoidalSystem *e, PetscReal lambda, PetscReal mu, PetscReal nu, PetscInt n, PetscInt p, PetscReal *val)
{
  PetscErrorCode ierr;
  PetscReal Enp, Inp;
  PetscInt signm;
  PetscInt signn;
  PetscFunctionBegin;

  if(mu < 0)
    signm = -1;
  else signm = 1;
  if(nu < 0)
    signn = -1;
  else signn = 1;


  ierr = CalcSolidInteriorHarmonic(e, lambda, mu, nu, n, p, &Enp);CHKERRQ(ierr);
  Inp = calcI(e, n, p, lambda, signm, signn);

  *val = (2*n + 1) * Enp * Inp;

  PetscLogFlops(4);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CalcSolidExteriorHarmonicVec"
PetscErrorCode CalcSolidExteriorHarmonicVec(EllipsoidalSystem* e, PetscInt nPoints, Vec ellPoints, PetscInt n, PetscInt p, Vec values)
{
  PetscErrorCode ierr;
  PetscReal Enp, Inp, Fnp;
  PetscReal lambda, mu, nu;
  PetscInt signm;
  PetscInt signn;
  const PetscScalar* ellPointsArray;
  PetscFunctionBegin;


  ierr = VecGetArrayRead(ellPoints, &ellPointsArray);CHKERRQ(ierr);
  for(PetscInt k=0; k<nPoints; ++k) {
    lambda = ellPointsArray[3*k+0];
    mu = ellPointsArray[3*k+1];
    nu = ellPointsArray[3*k+2];
    
    if(mu < 0)
      signm = -1;
    else signm = 1;
    if(nu < 0)
      signn = -1;
    else signn = 1;
    
    ierr = CalcSolidInteriorHarmonic(e, lambda, mu, nu, n, p, &Enp);CHKERRQ(ierr);
    Inp = calcI(e, n, p, lambda, signm, signn);
    Fnp = (2*n + 1) * Enp * Inp;
    ierr = VecSetValues(values, 1, &k, &Fnp, INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = VecAssemblyBegin(values);CHKERRQ(ierr); ierr = VecAssemblyEnd(values);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(ellPoints, &ellPointsArray);CHKERRQ(ierr);

  ierr = PetscLogFlops(4*nPoints);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "CalcCoulombEllCoefs"
/*
  GnpVals should be uninitialized. I
*/
PetscErrorCode CalcCoulombEllCoefs(EllipsoidalSystem* e, PetscInt nSource, Vec srcPoints, Vec srcCharges, PetscInt Nmax, Vec* GnpVals)
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

  
  ierr = VecCreateSeq(PETSC_COMM_SELF, nSource, &EnpVals);CHKERRQ(ierr);

  count = 0;
  for(PetscInt n=0; n<=Nmax; ++n) {
    for(PetscInt p=0; p<2*n+1; ++p)
      count++;
  }
  //ierr = VecDestroy(GnpVals);CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, count, GnpVals);CHKERRQ(ierr);
  
  ierr = VecGetArrayRead(srcCharges, &srcChargesArray);CHKERRQ(ierr);
  count = 0;
  for(PetscInt n=0; n<=Nmax; ++n) {
    for(PetscInt p=0; p<2*n+1; ++p) {
      normConstant = calcNormalization(e, n, p);
      
      ierr = CalcSolidInteriorHarmonicVec(e, nSource, srcPoints, n, p, EnpVals);
      
      ierr = VecGetArrayRead(srcPoints, &srcPointsArray);CHKERRQ(ierr);
      ierr = VecGetArrayRead(EnpVals, &EnpValsArray);CHKERRQ(ierr);
      Gnp = 0;
      for(PetscInt k=0; k<nSource; ++k) {
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

  ierr = PetscLogFlops(count * (6 + 3*nSource));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "CalcReactAndExtCoefsFromCoulomb"
PetscErrorCode CalcReactAndExtCoefsFromCoulomb(EllipsoidalSystem* e, PetscReal eps1, PetscReal eps2, PetscInt Nmax, Vec coulCoefs, Vec reacCoefs, Vec extCoefs)
{
  PetscErrorCode ierr;
  PetscReal Ea   , Ia   , Fa;
  PetscReal EaDer, IaDer, FaDer;
  PetscReal Gnp, Bnp, Cnp;
  PetscReal temp;
  const PetscScalar* coulCoefsArray;
  PetscInt count;
  PetscFunctionBegin;


  ierr = VecGetArrayRead(coulCoefs, &coulCoefsArray);CHKERRQ(ierr);
  count = 0;
  for(PetscInt n=0; n<=Nmax; ++n) {
    for(PetscInt p=0; p<2*n+1; ++p) {

      Gnp = coulCoefsArray[count];
      
      Ea    = calcLame(e, n, p, e->a, 1, 1);
      Ia    = calcI   (e, n, p, e->a, 1, 1);
      Fa    = (2*n + 1) * Ea    * Ia;
      EaDer = calcLameDerivative(e, n, p, e->a, 1, 1);
      IaDer = calcIDerivative   (e, n, p, e->a, 1, 1);
      FaDer = (2*n + 1) * (Ea*IaDer + EaDer*Ia);

      temp = (Fa/Ea)*(eps1 - eps2)/(eps1*eps2);
      temp /= (1 - (eps1/eps2)*((EaDer*Fa)/(FaDer*Ea)));
      Bnp  = temp*Gnp;

      Cnp = (eps1/eps2)*(EaDer/FaDer)*Bnp;
      Cnp += (Gnp/eps2);
      
      ierr = VecSetValues(reacCoefs, 1, &count, &Bnp, INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecSetValues(extCoefs, 1, &count, &Cnp, INSERT_VALUES);CHKERRQ(ierr);

      count++;
    }
  }

  ierr = VecAssemblyBegin(reacCoefs);CHKERRQ(ierr);
  ierr = VecAssemblyEnd  (reacCoefs);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(extCoefs);CHKERRQ(ierr);
  ierr = VecAssemblyEnd  (extCoefs);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(coulCoefs, &coulCoefsArray);CHKERRQ(ierr);


  ierr = PetscLogFlops(29*count);
  PetscFunctionReturn(0);
}


