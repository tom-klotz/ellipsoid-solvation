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
#include "../tanhsinh.h"


/*
EX3 - PLOTS CONVERGENCE OF A SINGLE NORMALIZATION CONSTANT
normConv.py PLOTS RESULT
*/


#undef __FUNCT__
#define __FUNCT__ "NormConstantIntFixed"
PetscErrorCode NormConstantIntFixed(EllipsoidalSystem *e, PetscInt n, PetscInt p, PetscInt prec, PetscInt nPts, PetscReal *intVals, PetscReal *normConst)
{
  PetscErrorCode ierr;
  PetscInt flopCount;
  PetscFunctionBegin;
  flopCount = 0;

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
  intVals[0] = integrals[0];
  ierr = DEQuad((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, e->hp_h, e->hp_k, prec, nPts, integrals+1, &ctx2);
  intVals[1] = integrals[1];
  ierr = DEQuad((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, mpfrzero, e->hp_h, prec, nPts, integrals+2, &ctx3);
  intVals[2] = integrals[2];
  ierr = DEQuad((PetscErrorCode (*)(mpfr_t*, mpfr_t*, void*)) normFunction1, mpfrzero, e->hp_h, prec, nPts, integrals+3, &ctx4);
  intVals[3] = integrals[3];
  *normConst = 8.0*(integrals[2]*integrals[1] - integrals[0]*integrals[3]); flopCount += 4;

  mpfr_clears(mpfrzero, mpfrone, NULL);
  
  ierr = PetscLogFlops(flopCount);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT_ "normConvergence"
PetscErrorCode normConvergence(PetscInt n, PetscInt p)
{
  PetscErrorCode ierr;
  
  const PetscInt MPFR_PREC = 64;
  const PetscInt NUM_SOLUTIONS = 80;
  const PetscInt POINTS_MIN = 4;
  const PetscInt POINTS_STEP = 2;
  const PetscInt prec = 16;
  const PetscReal a = 3.0;
  const PetscReal b = 2.0;
  const PetscReal c = 1.0;
  mpfr_t solExact;
  mpfr_t intExact[4];
  PetscReal solutions[NUM_SOLUTIONS];
  mpfr_t errors   [NUM_SOLUTIONS];
  PetscReal intSols[4*NUM_SOLUTIONS];
  mpfr_t intErrors[4*NUM_SOLUTIONS];
  PetscReal sol2, sol3, sol4;
  PetscLogEvent flopCounts[NUM_SOLUTIONS];
  EllipsoidalSystem e;
  PetscInt i;
  PetscEventPerfInfo info;
  //PetscInt n = 7;
  //PeatscInt p = 4;
  
  PetscFunctionBegin;


  mpfr_set_default_prec(4*MPFR_PREC);
  initEllipsoidalSystem(&e, a, b, c, MPFR_PREC);


  mpfr_init(solExact);
  for(PetscInt k=0; k < 4; ++k)
    mpfr_init(intExact[k]);
  for(PetscInt k=0; k < NUM_SOLUTIONS; ++k )
    mpfr_init(errors[k]);
  for(PetscInt k=0; k < 4*NUM_SOLUTIONS; ++k)
    mpfr_init(intErrors[k]);
  
  
  /* calculate "exact" solution */
  ierr = calcNormalization2MPFR(&e, n, p, intExact, &solExact);CHKERRQ(ierr);
  printf("old norm constant: %15.15f\n", solExact);

  /* calculate approximate solutions and record flops */
  char text[40] = "%d points";
  char sText[40];
  for(i=0; i<NUM_SOLUTIONS; ++i) {
    sprintf(sText, text, POINTS_MIN + POINTS_STEP*i);
    ierr = PetscLogEventRegister(sText, 0, flopCounts+i);CHKERRQ(ierr);
    ierr = PetscLogEventBegin(flopCounts[i], 0, 0, 0, 0);CHKERRQ(ierr);
    ierr = NormConstantIntFixed(&e, n, p, prec, POINTS_MIN + POINTS_STEP*i, intSols+4*i, solutions+i);CHKERRQ(ierr);
    ierr = PetscLogEventEnd(flopCounts[i], 0, 0, 0, 0);CHKERRQ(ierr);
    printf("new norm constant: %15.15f\n", solutions[i]);
  }
  /* calculate errors */
  for(i=0; i<NUM_SOLUTIONS; ++i) {
    //errors[i] = PetscAbsReal((solExact - solutions[i])/solExact);
    mpfr_sub_d(errors[i], solExact, solutions[i], MPFR_RNDN);
    mpfr_div(errors[i], errors[i], solExact, MPFR_RNDN);
    mpfr_abs(errors[i], errors[i], MPFR_RNDN);
    //intErrors[4*i+0] = PetscAbsReal((intExact[0] - intSols[4*i+0])/intExact[0]);
    mpfr_sub_d(intErrors[4*i+0], intExact[0], intSols[4*i+0], MPFR_RNDN);
    mpfr_div(intErrors[4*i+0], intErrors[4*i+0], intExact[0], MPFR_RNDN);
    mpfr_abs(intErrors[4*i+0], intErrors[4*i+0], MPFR_RNDN);
    //intErrors[4*i+1] = PetscAbsReal((intExact[1] - intSols[4*i+1])/intExact[1]);
    mpfr_sub_d(intErrors[4*i+1], intExact[1], intSols[4*i+1], MPFR_RNDN);
    mpfr_div(intErrors[4*i+1], intErrors[4*i+1], intExact[1], MPFR_RNDN);
    mpfr_abs(intErrors[4*i+1], intErrors[4*i+1], MPFR_RNDN);
    //intErrors[4*i+2] = PetscAbsReal((intExact[2] - intSols[4*i+2])/intExact[2]);
    mpfr_sub_d(intErrors[4*i+2], intExact[2], intSols[4*i+2], MPFR_RNDN);
    mpfr_div(intErrors[4*i+2], intErrors[4*i+2], intExact[2], MPFR_RNDN);
    mpfr_abs(intErrors[4*i+2], intErrors[4*i+2], MPFR_RNDN);
    //intErrors[4*i+3] = PetscAbsReal((intExact[3] - intSols[4*i+3])/intExact[3]);
    mpfr_sub_d(intErrors[4*i+3], intExact[3], intSols[4*i+3], MPFR_RNDN);
    mpfr_div(intErrors[4*i+3], intErrors[4*i+3], intExact[3], MPFR_RNDN);
    mpfr_abs(intErrors[4*i+3], intErrors[4*i+3], MPFR_RNDN);
    printf("errors[%d] = %15.15f\n", i, errors[i]);
  }
  
  FILE *fp1 = fopen("out/normInt1Prec.txt", "w");
  FILE *fp2 = fopen("out/normInt2Prec.txt", "w");
  FILE *fp3 = fopen("out/normInt3Prec.txt", "w");
  FILE *fp4 = fopen("out/normInt4Prec.txt", "w");
  FILE *fp = fopen("out/normWorkPrec.txt", "w");
  fprintf(fp, "points flops error\n");
  fprintf(fp1, "points error\n");
  fprintf(fp2, "points error\n");
  fprintf(fp3, "points error\n");
  fprintf(fp4, "points error\n");
  for(i=0; i<NUM_SOLUTIONS; ++i) {
    ierr = PetscLogEventGetPerfInfo(PETSC_DETERMINE, flopCounts[i], &info);CHKERRQ(ierr);
    fprintf(fp, "%d %4.4e %4.4e\n", POINTS_MIN + POINTS_STEP*i, info.flops, mpfr_get_d(errors[i], MPFR_RNDN));
    fprintf(fp1, "%d %4.4e\n", POINTS_MIN + POINTS_STEP*i, mpfr_get_d(intErrors[4*i+0], MPFR_RNDN));
    fprintf(fp2, "%d %4.4e\n", POINTS_MIN + POINTS_STEP*i, mpfr_get_d(intErrors[4*i+1], MPFR_RNDN));
    fprintf(fp3, "%d %4.4e\n", POINTS_MIN + POINTS_STEP*i, mpfr_get_d(intErrors[4*i+2], MPFR_RNDN));
    fprintf(fp4, "%d %4.4e\n", POINTS_MIN + POINTS_STEP*i, mpfr_get_d(intErrors[4*i+3], MPFR_RNDN));
  }
  fclose(fp);
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
  mpfr_clear(solExact);

  PetscFunctionReturn(0);

}

#undef __FUNCT__
#define __FUNCT__ "main"
PetscErrorCode main(int argc, char **argv)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;

  ierr = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRQ(ierr);
  ierr = PetscLogDefaultBegin();CHKERRQ(ierr);

  ierr = normConvergence(21,22);CHKERRQ(ierr);
  
  ierr = PetscFinalize();

  PetscFunctionReturn(0);
}
