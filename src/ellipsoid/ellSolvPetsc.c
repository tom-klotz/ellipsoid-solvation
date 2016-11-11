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
#define __FUNCT__ "InitEllipsoidalAndConvertPoints"
PetscErrorCode InitEllipsoidalAndConvertPoints(EllipsoidalSystem *e, PetscReal a, PetscReal b, PetscReal c, Vec xyz, Vec ell)
{
  PetscErrorCode ierr;
  PetscLogEvent initEvent;
  PetscFunctionBegin;

  ierr = PetscLogEventRegister("Init ellipsoidal system", 0, &initEvent);CHKERRQ(ierr);
  ierr = PetscLogEventBegin(initEvent, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = initEllipsoidalSystem(e, a, b, c);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(initEvent, 0, 0, 0, 0);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "HowMany"
PetscErrorCode HowMany(PetscInt N, PetscInt *num)
{
  PetscErrorCode ierr;
  PetscInt n, p;
  PetscFunctionBegin;
  *num = 0;
  for(n = 0; n <= N; ++n) {
    for(p=0; p < 2*n + 1; ++p)
      *num += 1;
  }
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "CalcEllipsoidFreeEnergy"
PetscErrorCode CalcEllipsoidFreeEnergy(EllipsoidalSystem *e, PetscReal eps1, PetscReal eps2, PetscInt nSrc, Vec srcXYZ, Vec srcMag, PetscReal tol, PetscInt Nmax, Vec tarSol, PetscReal *freeE)
{
  PetscErrorCode ierr;
  PetscInt flopCount;
  //EllipsoidalSystem e;
  PetscReal x, y, z;
  PetscInt n, p, npsize;
  Vec srcEll, tarEll;
  PetscScalar *tarEllArray;
  Vec coulCoefs, reactCoefs, extCoefs;
  const PetscScalar *coulCoefsArray, *reactCoefsArray, *extCoefsArray;
  const PetscScalar *srcMagArray;
  Vec EnpVals, FnpVals;
  const PetscScalar *EnpValsArray, *FnpValsArray;
  PetscScalar *tarSolArray;
  PetscFunctionBegin;
  flopCount = 0;
  
  /* init ellipsoidal system */
  //ierr = initEllipsoidalSystem(&e, a, b, c);CHKERRQ(ierr);
  
  /* initialize expansion vectors */
  ierr = HowMany(Nmax, &npsize);CHKERRQ(ierr);
  //ierr = VecCreateSeq(PETSC_COMM_SELF, npsize, &coulCoefs);CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, npsize, &reactCoefs);CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, npsize, &extCoefs);CHKERRQ(ierr);
  
  /* initialize point vectors */
  ierr = VecDuplicate(srcXYZ, &srcEll);CHKERRQ(ierr);
  ierr = VecDuplicate(srcEll, &tarEll);CHKERRQ(ierr);
  
  /* convert points from cartesian to ellipsoidal */
  ierr = CartesianToEllipsoidalVec(e, srcXYZ, srcEll);CHKERRQ(ierr);
  ierr = VecCopy(srcEll, tarEll);CHKERRQ(ierr); //source+target same

  /* Calculate Gnp */
  ierr = CalcCoulombEllCoefs(e, nSrc, srcEll, srcMag, Nmax, &coulCoefs);CHKERRQ(ierr);
  /* calculate Bnp and Cnp */
  ierr = CalcReactAndExtCoefsFromCoulomb(e, eps1, eps2, Nmax, coulCoefs, reactCoefs, extCoefs);

  // create EnpVals and FnpVals vectors
  ierr = VecCreateSeq(PETSC_COMM_SELF, nSrc, &EnpVals);CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, nSrc, &FnpVals);CHKERRQ(ierr);
  
  // get read-only pointers for expansion coefficients
  ierr = VecGetArrayRead(coulCoefs, &coulCoefsArray);CHKERRQ(ierr);
  ierr = VecGetArrayRead(reactCoefs, &reactCoefsArray);CHKERRQ(ierr);
  ierr = VecGetArrayRead(extCoefs, &extCoefsArray);CHKERRQ(ierr);
  // get write pointer for solution vector
  ierr = VecZeroEntries(tarSol);CHKERRQ(ierr);
  ierr = VecGetArray(tarSol, &tarSolArray);CHKERRQ(ierr);
  PetscInt ind = 0;
  for(n=0; n <= Nmax; ++n) {
    for(p=0; p < 2*n+1; ++p) {
      PetscReal Gnp =  coulCoefsArray[ind];
      PetscReal Cnp =   extCoefsArray[ind];
      PetscReal Bnp = reactCoefsArray[ind];
      
      ierr = CalcSolidInteriorHarmonicVec(e, tarEll, n, p, EnpVals);CHKERRQ(ierr);
      //ierr = CalcSolidExteriorHarmonicVec(&e, tarEll, n, p, FnpVals);CHKERRQ(ierr);

      ierr = VecGetArrayRead(EnpVals, &EnpValsArray);CHKERRQ(ierr);
      //ierr = VecGetArrayRead(FnpVals, &FnpValsArray);CHKERRQ(ierr);
      for(PetscInt k=0; k < nSrc; ++k) {
	PetscReal Enp = EnpValsArray[k];
	//PetscReal Fnp = FnpValsArray[k];
	tarSolArray[k] += Bnp*Enp; flopCount += 2;
      }
      ierr = VecRestoreArrayRead(EnpVals, &EnpValsArray);CHKERRQ(ierr);
      //ierr = VecRestoreArrayRead(FnpVals, &FnpValsArray);CHKERRQ(ierr);

      ind++;
    }
  }
  // get read-only pointers for expansion coefficients
  ierr = VecRestoreArrayRead(coulCoefs, &coulCoefsArray);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(reactCoefs, &reactCoefsArray);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(extCoefs, &extCoefsArray);CHKERRQ(ierr);
  *freeE = 0;
  ierr = VecGetArrayRead(srcMag, &srcMagArray);CHKERRQ(ierr);
  /* calculate free energy */
  for(PetscInt k=0; k < nSrc; ++k) {
    *freeE += srcMagArray[k]*tarSolArray[k];
  }
  ierr = VecRestoreArrayRead(srcMag, &srcMagArray);CHKERRQ(ierr);
  ierr = VecRestoreArray(tarSol, &tarSolArray);CHKERRQ(ierr);
  ierr = PetscLogFlops(flopCount);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CalcEllipsoidTester"
PetscErrorCode CalcEllipsoidTester(PetscReal a, PetscReal b, PetscReal c, PetscReal eps1, PetscReal eps2, PetscInt nSrc, Vec srcXYZ, Vec srcMag, PetscInt nTar, Vec tarXYZ, PetscInt Nmax, Vec tarSol)
{
  PetscErrorCode ierr;
  EllipsoidalSystem e;
  PetscReal x, y, z;
  PetscInt n, p, npsize;
  Vec srcEll, tarEll;
  PetscScalar *tarEllArray;
  Vec coulCoefs, reactCoefs, extCoefs;
  const PetscScalar *coulCoefsArray, *reactCoefsArray, *extCoefsArray;
  Vec EnpVals, FnpVals;
  const PetscScalar *EnpValsArray, *FnpValsArray;
  PetscScalar *tarSolArray;
  PetscFunctionBegin;


  /* init ellipsoidal system */
  ierr = initEllipsoidalSystem(&e, a, b, c);CHKERRQ(ierr);

  /* initialize expansion vectors */
  ierr = HowMany(Nmax, &npsize);CHKERRQ(ierr);
  //ierr = VecCreateSeq(PETSC_COMM_SELF, npsize, &coulCoefs);CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, npsize, &reactCoefs);CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, npsize, &extCoefs);CHKERRQ(ierr);
  
  /* initialize ellipsoid vectors */
  ierr = VecDuplicate(srcXYZ, &srcEll);CHKERRQ(ierr);
  ierr = VecDuplicate(tarXYZ, &tarEll);CHKERRQ(ierr);
  
  /* convert points from cartesian to ellipsoidal */
  ierr = CartesianToEllipsoidalVec(&e, srcXYZ, srcEll);CHKERRQ(ierr);
  ierr = CartesianToEllipsoidalVec(&e, tarXYZ, tarEll);CHKERRQ(ierr);
  

  /* Calculate Gnp */
  ierr = CalcCoulombEllCoefs(&e, nSrc, srcEll, srcMag, Nmax, &coulCoefs);CHKERRQ(ierr);
  //printf("\n\n##################################\n########## COUL COEFS #############\n###########################################\n\n");
  //printf("nmax: %d\n", Nmax);
  //ierr = VecView(coulCoefs, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
  //printf("\n\n##################################\n########## OVER #############\n###########################################\n\n");
  /* calculate Bnp and Cnp */
  ierr = CalcReactAndExtCoefsFromCoulomb(&e, eps1, eps2, Nmax, coulCoefs, reactCoefs, extCoefs);

  // create EnpVals and FnpVals vectors
  ierr = VecCreateSeq(PETSC_COMM_SELF, nTar, &EnpVals);CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, nTar, &FnpVals);CHKERRQ(ierr);
  
  // get read-only pointers for expansion coefficients
  ierr = VecGetArrayRead(coulCoefs, &coulCoefsArray);CHKERRQ(ierr);
  ierr = VecGetArrayRead(reactCoefs, &reactCoefsArray);CHKERRQ(ierr);
  ierr = VecGetArrayRead(extCoefs, &extCoefsArray);CHKERRQ(ierr);
  // get write pointer for solution vector
  ierr = VecZeroEntries(tarSol);CHKERRQ(ierr);
  ierr = VecGetArray(tarSol, &tarSolArray);CHKERRQ(ierr);
  PetscInt ind = 0;
  for(n=0; n <= Nmax; ++n) {
    printf("n: %d\n", n);
    for(p=0; p < 2*n+1; ++p) {
      PetscReal Gnp =  coulCoefsArray[ind];
      PetscReal Cnp =   extCoefsArray[ind];
      PetscReal Bnp = reactCoefsArray[ind];
      
      ierr = CalcSolidInteriorHarmonicVec(&e, tarEll, n, p, EnpVals);CHKERRQ(ierr);
      ierr = CalcSolidExteriorHarmonicVec(&e, tarEll, n, p, FnpVals);CHKERRQ(ierr);

      ierr = VecGetArrayRead(EnpVals, &EnpValsArray);CHKERRQ(ierr);
      ierr = VecGetArrayRead(FnpVals, &FnpValsArray);CHKERRQ(ierr);
      for(PetscInt k=0; k < nTar; ++k) {
	if(n==0 && p==0) {
	  tarSolArray[k] = 0;
	}
	PetscReal lambda;
	PetscInt index = 3*k+0;
	ierr = VecGetValues(tarEll, 1, &index, &lambda);CHKERRQ(ierr);
	
	PetscReal Enp = EnpValsArray[k];
	PetscReal Fnp = FnpValsArray[k];
	if(PetscAbsReal(lambda) <= a)
	  tarSolArray[k] += Bnp*Enp; //(Gnp/eps1)*Fnp;
	else
	  tarSolArray[k] += Cnp*Fnp;//(Gnp/eps1)*Fnp;
	
      }
      ierr = VecRestoreArrayRead(EnpVals, &EnpValsArray);CHKERRQ(ierr);
      ierr = VecRestoreArrayRead(FnpVals, &FnpValsArray);CHKERRQ(ierr);

      ind++;
    }
  }
  
  PetscFunctionReturn(0);
}

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

  //PetscLogStage stageIntSolids, stageExtSolids;
  PetscLogEvent *INTERIOR_FLOPS;
  PetscLogEvent *EXTERIOR_FLOPS;
  //PetscLogDouble user_log_flops;
  PetscFunctionBegin;

  
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
  PetscInt howmany = 0;
  for(n=0; n <= Nmax; ++n) {
    for(p=0; p < 2*n+1; ++p)
      howmany++;
  }

  ierr = PetscMalloc1(sizeof(PetscLogEvent)*howmany, &INTERIOR_FLOPS);CHKERRQ(ierr);
  ierr = PetscMalloc1(sizeof(PetscLogEvent)*howmany, &EXTERIOR_FLOPS);CHKERRQ(ierr);
  char intStr[40] = "Interior Harmonic %d";
  char extStr[40] = "Exterior Harmonic %d";
  char indvIntStr[40];
  char indvExtStr[40];
  for(n=0; n<howmany; ++n) {
    sprintf(indvIntStr, intStr, n);
    sprintf(indvExtStr, extStr, n);
    ierr = PetscLogEventRegister(indvIntStr, 0, INTERIOR_FLOPS+n);CHKERRQ(ierr);
    ierr = PetscLogEventRegister(indvExtStr, 0, EXTERIOR_FLOPS+n);CHKERRQ(ierr);
  }
  PetscInt ind = 0;
  PetscInt intInd;
  PetscInt extInd;
  for(n=0; n <= Nmax; ++n) {
    printf("%d/%d\n", n, Nmax);
    for(p=0; p < 2*n+1; ++p) {
      
      PetscReal Gnp =  coulCoefsArray[ind];
      PetscReal Cnp =   extCoefsArray[ind];
      PetscReal Bnp = reactCoefsArray[ind];
      // calc solid harmonics for interior,exterior points


      //ierr = PetscLogStagePush(stageIntSolids);CHKERRQ(ierr);
      ierr = PetscLogEventBegin(INTERIOR_FLOPS[ind],0,0,0,0);CHKERRQ(ierr);
      ierr = CalcSolidInteriorHarmonicVec(&e, tarIntEll, n, p, EnpVals);CHKERRQ(ierr);
      ierr = PetscLogEventEnd  (INTERIOR_FLOPS[ind],0,0,0,0);CHKERRQ(ierr);
      //ierr = PetscLogStagePop();
      
      //ierr = PetscLogStagePush(stageExtSolids);CHKERRQ(ierr);
      ierr = PetscLogEventBegin(EXTERIOR_FLOPS[ind],0,0,0,0);CHKERRQ(ierr);
      ierr = CalcSolidExteriorHarmonicVec(&e, tarExtEll, n, p, FnpVals);CHKERRQ(ierr);
      ierr = PetscLogEventEnd(EXTERIOR_FLOPS[ind],0,0,0,0);CHKERRQ(ierr);
      //ierr = PetscLogStagePop();CHKERRQ(ierr);
      
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
	  //targetSolArray[k] += 1;
	  intInd++;
	}
	else { // exterior points
	  //printf("F%d%d: %15.15f\nCnp %d%d: %15.15f\n\n", n, p, n, p, FnpValsArray[intInd], extCoefsArray[ind]);
	  //targetSolArray[k] += extCoefsArray[ind]*FnpValsArray[extInd];
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
  ierr = VecRestoreArrayRead(coulCoefs , &coulCoefsArray );CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(reactCoefs, &reactCoefsArray);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(extCoefs  , &extCoefsArray  );CHKERRQ(ierr);


  //output flops to file
  FILE *fpext = fopen("flopsext.txt", "w");
  FILE *fpint = fopen("flopsint.txt", "w");
  PetscStageLog stageLog;  
  ierr = PetscLogGetStageLog(&stageLog);CHKERRQ(ierr);

  /*
  for(int i=0; i<howmany; ++i) {
    fprintf(fpext, "Ext %d = %.4e\n", i, stageLog->stageInfo[stageExtSolids].eventLog->eventInfo[EXTERIOR_FLOPS[i]].flops);
    fprintf(fpint, "Int %d = %.4e\n", i, stageLog->stageInfo[stageIntSolids].eventLog->eventInfo[INTERIOR_FLOPS[i]].flops);
  }
  */


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

  calcLame(e, n, p, lambda, signm, signn, &eL);
  calcLame(e, n, p, mu, signm, signn, &eM);
  calcLame(e, n, p, nu, signm, signn, &eN);
  
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
PetscErrorCode CalcSolidInteriorHarmonicVec(EllipsoidalSystem* e, Vec ellPoints, PetscInt n, PetscInt p, Vec values)
{
  PetscErrorCode ierr;
  PetscInt signm = 1;
  PetscInt signn = 1;
  PetscInt nPoints;
  PetscReal lambda, mu, nu;
  PetscReal eL, eM, eN;
  PetscReal Enp;
  const PetscScalar* ellPointsArray;
  PetscFunctionBegin;
  ierr = VecGetSize(ellPoints, &nPoints);CHKERRQ(ierr);
  nPoints = nPoints/3;
  
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
    ierr = calcLame(e, n, p, lambda, signm, signn, &eL);
    ierr = calcLame(e, n, p, mu, signm, signn, &eM);
    ierr = calcLame(e, n, p, nu, signm, signn, &eN);
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
  ierr = calcI(e, n, p, lambda, signm, signn, &Inp);CHKERRQ(ierr);

  *val = (2*n + 1) * Enp * Inp;

  PetscLogFlops(4);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CalcSolidExteriorHarmonicVec"
PetscErrorCode CalcSolidExteriorHarmonicVec(EllipsoidalSystem* e, Vec ellPoints, PetscInt n, PetscInt p, Vec values)
{
  PetscErrorCode ierr;
  PetscReal Enp, Inp, Fnp;
  PetscReal lambda, mu, nu;
  PetscInt signm;
  PetscInt signn;
  PetscInt nPoints;
  const PetscScalar* ellPointsArray;
  PetscFunctionBegin;
  ierr = VecGetSize(ellPoints, &nPoints);CHKERRQ(ierr);
  nPoints = nPoints/3;

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
    ierr = calcI(e, n, p, lambda, signm, signn, &Inp);CHKERRQ(ierr);
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
      ierr = calcNormalization(e, n, p, &normConstant);CHKERRQ(ierr);
      
      ierr = CalcSolidInteriorHarmonicVec(e, srcPoints, n, p, EnpVals);
      
      ierr = VecGetArrayRead(EnpVals, &EnpValsArray);CHKERRQ(ierr);
      Gnp = 0;
      for(PetscInt k=0; k<nSource; ++k) {
	Gnp += srcChargesArray[k]*EnpValsArray[k];
      }
      ierr = VecRestoreArrayRead(EnpVals, &EnpValsArray);CHKERRQ(ierr);
      Gnp *= (4*PETSC_PI)/((2.0*n+1.0)*normConstant);
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
      
      ierr = calcLame(e, n, p, e->a, 1, 1, &Ea);CHKERRQ(ierr);
      ierr = calcI   (e, n, p, e->a, 1, 1, &Ia);CHKERRQ(ierr);
      Fa    = (2*n + 1) * Ea    * Ia;
      calcLameDerivative(e, n, p, e->a, 1, 1, &EaDer);
      calcIDerivative   (e, n, p, e->a, 1, 1, &IaDer);
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


