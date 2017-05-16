#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpfr.h>
#include <petsc.h>

#include "../ellipsoid/ellipsoid.h"
#include "../ellipsoid/ellSolv.h"
#include "../sphere/sphere.h"
#include "../sphere/sphSolv.h"
#include "../constants.h"
#include "testFunctions.h"

//compares ellipsoidal and spherical
//has ellipsoidal distribution of charges just inside ellipsoidal boundary
//computes the solution at two point, one inside smallest enclosing sphere where spherical diverges
//other solution point is outside of sphere and ellipse to compare convergence
#undef __FUNCT__
#define __FUNCT__ "main"
PetscErrorCode main(int argc, char **argv)
{
  PetscErrorCode ierr;
  //semi-axes of ellipsoid
  PetscReal a = 3.0;
  PetscReal b = 2.5;
  PetscReal c = 1.2;
  //sphere radius
  PetscReal sphRad = 3;
  //dielecric permitivities both set to 1
  const PetscReal eps1 = 1.0;
  const PetscReal eps2 = 1.0;
  //charges placed on ellipsoid with sem-axes smaller by epsilon
  PetscReal epsilon = .5;
  //interior charge XYZ and magnitude
  const PetscInt SQRT_NUM_CHARGES = 3;
  Vec chargeXYZ;
  Vec chargeMag;
  //maximum expansion order for test
  const PetscInt MAX_N = 25;
  //solutions vector
  const PetscInt NUM_SOL_PTS = 2;
  Vec ellSolutions[MAX_N], sphSolutions[MAX_N];
  Vec ellError[MAX_N], sphError[MAX_N];
  Vec exactSolution;
  Vec solXYZ;
  //solution point inside of sphere
  PetscReal sol1X = .1;
  PetscReal sol1Y = .1;
  PetscReal sol1Z = 1.5;
  //solution point outside of sphere
  PetscReal sol2X = 2.1;
  PetscReal sol2Y = 2.6;
  PetscReal sol2Z = 1;
  //pointer for accessing petsc vectors
  PetscScalar *vecPt, *vecPt2;
  PetscFunctionBeginUser;

  //initialize petsc
  ierr = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRQ(ierr);
  ierr = PetscLogDefaultBegin();CHKERRQ(ierr);

  //create, set charge magnitudes to 1
  ierr = VecCreateSeq(PETSC_COMM_SELF, SQRT_NUM_CHARGES*SQRT_NUM_CHARGES, &chargeMag);CHKERRQ(ierr);
  ierr = VecSet(chargeMag, 1.0);CHKERRQ(ierr);

  //create, set charge XYZ values to be on a-eps,b-eps,c-eps ellipsoid
  ierr = VecCreateSeq(PETSC_COMM_SELF, 3*SQRT_NUM_CHARGES*SQRT_NUM_CHARGES, &chargeXYZ);CHKERRQ(ierr);
  ierr = VecGetArray(chargeXYZ, &vecPt);CHKERRQ(ierr);
  PetscReal pu, pv;
  for(PetscInt i=0; i<SQRT_NUM_CHARGES; ++i) {
    for(PetscInt j=0; j<SQRT_NUM_CHARGES; ++j) {
      pu = i*2.0*PETSC_PI/(SQRT_NUM_CHARGES-1);
      pv = j*PETSC_PI/(SQRT_NUM_CHARGES-1);
      PetscInt ind = i*SQRT_NUM_CHARGES + j;
      //there is a slight issue with on-axes coordinate transforms which is what the +.01 is for
      vecPt[3*ind + 0] = (a-epsilon)*PetscCosReal(pu)*PetscSinReal(pv) +.01;
      vecPt[3*ind + 1] = (b-epsilon)*PetscSinReal(pu)*PetscSinReal(pv) + .01;
      vecPt[3*ind + 2] = (c-epsilon)*PetscCosReal(pv) + .01;
      printf("pu: %3.3f pv: %3.3f\n", pu, pv);
    }
  }
  ierr = VecRestoreArray(chargeXYZ, &vecPt);CHKERRQ(ierr);

  //set solution points vector
  ierr = VecCreateSeq(PETSC_COMM_SELF, 3*NUM_SOL_PTS, &solXYZ);CHKERRQ(ierr);
  ierr = VecGetArray(solXYZ, &vecPt);CHKERRQ(ierr);
  vecPt[0] = sol1X; vecPt[1] = sol1Y; vecPt[2] = sol1Z;
  vecPt[3] = sol2X; vecPt[4] = sol2Y; vecPt[5] = sol2Z;
  ierr = VecRestoreArray(solXYZ, &vecPt);CHKERRQ(ierr);
  
  //caulculate exact solution
  ierr = VecCreateSeq(PETSC_COMM_SELF, NUM_SOL_PTS, &exactSolution);CHKERRQ(ierr);
  ierr = CoulombExact(eps1, chargeXYZ, chargeMag, solXYZ, exactSolution);CHKERRQ(ierr);

  //initialize solutions and error vectors
  for(PetscInt k=0; k<MAX_N; ++k) {
    ierr = VecCreateSeq(PETSC_COMM_SELF, NUM_SOL_PTS, ellSolutions+k);CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF, NUM_SOL_PTS, sphSolutions+k);CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF, NUM_SOL_PTS, ellError+k);CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF, NUM_SOL_PTS, sphError+k);CHKERRQ(ierr);
  }

  //compute solutions for expansions up to order MAX_N
  for(PetscInt k=1; k<=MAX_N; ++k) {
    printf("N = %d\n", k);
    printf("calculating ellipsoidal solution...\n");
    ierr = CalcEllipsoidCoulombPotential(a, b, c, eps1, eps2, SQRT_NUM_CHARGES*SQRT_NUM_CHARGES, chargeXYZ, chargeMag, NUM_SOL_PTS, solXYZ, k, ellSolutions[k-1]);CHKERRQ(ierr);
    printf("calculating spherical solution...\n");
    ierr = CalcSphericalCoulombPotential(1.0, eps1, eps2, SQRT_NUM_CHARGES*SQRT_NUM_CHARGES, chargeXYZ, chargeMag, NUM_SOL_PTS, solXYZ, k, sphSolutions[k-1]);CHKERRQ(ierr);
  }

  //compute errors
  printf("Computing errors...\n");
  for(PetscInt k=0; k<MAX_N; ++k) {
    ierr = VecCopy(exactSolution, ellError[k]);CHKERRQ(ierr);
    ierr = VecCopy(exactSolution, sphError[k]);CHKERRQ(ierr);
    ierr = VecAXPY(ellError[k], -1.0, ellSolutions[k]);CHKERRQ(ierr);
    ierr = VecAXPY(sphError[k], -1.0, sphSolutions[k]);CHKERRQ(ierr);
    ierr = VecAbs(ellError[k]);CHKERRQ(ierr);
    ierr = VecAbs(sphError[k]);CHKERRQ(ierr);
  }

  FILE *fp = fopen("out/ex1.txt", "w");
  
  
  for(PetscInt k=0; k<MAX_N; ++k) {
    ierr = VecGetArrayRead(ellError[k], &vecPt);
    ierr = VecGetArrayRead(sphError[k], &vecPt2);
    printf("Error[%d] = %4.4e\n", k, vecPt[0]);
    fprintf(fp, "%4.4e %4.4e %4.4e %4.4e\n", vecPt[0], vecPt2[0], vecPt[1], vecPt2[1]);
    ierr = VecRestoreArrayRead(ellError[k], &vecPt);
    ierr = VecRestoreArrayRead(sphError[k], &vecPt2);
  }
  printf("Exact: \n");
  ierr = VecView(exactSolution, PETSC_VIEWER_STDOUT_SELF);
  printf("Ell: \n");
  ierr = VecView(ellSolutions[MAX_N-1], PETSC_VIEWER_STDOUT_SELF);
  
  PetscFunctionReturn(0);
}
